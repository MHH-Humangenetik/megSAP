<?php
/** 
	@page data_setup
*/

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

//parse command line arguments
$parser = new ToolBase("data_setup", "Creates a local copy of the reference genome, annotation data and apptainer container. They are heavily used during data analysis and should not be accessed via the network.");
$parser->addString("build", "Genome build.", true, "GRCh38");
$parser->addFlag("check", "Check annotation files.");
extract($parser->parse($argv));

//make sure there are no singularity left-overs
if (file_exists("/etc/singularity/"))
{
	trigger_error("The server contains singularity left-overs: cannot execute megSAP until you delete /etc/singularity/!", E_USER_ERROR);
}

//init reference genome and annotation data
$data_folder = get_path("data_folder");
$local_data = get_path("local_data");
$rsync  = "rsync --recursive --no-perms --no-acls --omit-dir-times --no-group --no-owner --chmod=ugo=rwX --copy-links --size-only";

//determine DB files
$db_files = array("/dbs/CADD/CADD_SNVs_1.7_GRCh38.vcf.gz", "/dbs/CADD/CADD_InDels_1.7_GRCh38.vcf.gz", "/dbs/REVEL/REVEL_1.3.vcf.gz", "/dbs/AlphaMissense/AlphaMissense_hg38.vcf.gz", "/dbs/gnomAD/gnomAD_genome_v4.1_GRCh38.vcf.gz", "/dbs/gnomAD/gnomAD_genome_v3.1.mito_GRCh38.vcf.gz", "/dbs/RepeatMasker/RepeatMasker_GRCh38.bed", "/dbs/ClinVar/clinvar_20250128_converted_GRCh38.vcf.gz", "/dbs/phyloP/hg38.phyloP100way.bw", "/dbs/SpliceAI/spliceai_scores_2024_08_26_GRCh38.vcf.gz");
$omim =  "/dbs/OMIM/omim.bed";
if (file_exists($data_folder.$omim)) $db_files[] = $omim; //optional
$hgmd =  "/dbs/HGMD/HGMD_PRO_2024_4_fixed.vcf.gz";
if (file_exists($data_folder.$hgmd)) $db_files[] = $hgmd; //optional

######################### reference genome #########################
$genome_folder = "{$data_folder}/genomes/";
print "### Reference genome ###\n";
print "from: {$genome_folder}\n";
print "to  : {$local_data}\n";
print "\n";

//make local folder
if (!file_exists($local_data))
{
	if (!mkdir($local_data))
	{
		trigger_error("Could not create local data folder '{$local_data}'!", E_USER_ERROR);
	}
	if (!chmod($local_data, 0777))
	{
		trigger_error("Could not change privileges of local data folder '{$local_data}'!", E_USER_ERROR);
	}
}

//remove outdated genome FASTA files
if (file_exists("{$local_data}/{$build}.fa.md5"))
{
	exec("diff {$local_data}/{$build}.fa.md5 {$genome_folder}/{$build}.fa.md5", $output, $code);
	if ($code==0)
	{
		print "MD5 checksums of genome '{$build}.fa.md5' match.\n";
	}
	else
	{
		print "MD5 checksums of genome '{$build}.fa.md5' differ. Deleting old data!\n";
		exec2("rm -rf {$local_data}/{$build}.dict {$local_data}/{$build}.fa*");
	}
}

//wait sync is done by another processes
$pid_file = "{$local_data}/megSAP_data_setup_{$build}.txt";
print "PID: $pid_file\n";
if (file_exists($pid_file))
{
	$pid_old = (int)trim(file_get_contents($pid_file));
	$iter = 0;
	while (posix_getpgid($pid_old)!==FALSE)
	{
		++$iter;
		if ($iter>180) break; //wait for 3h max
		print "Process with PID {$pid_old} is already syncing. Waiting one minute...\n";
		sleep(60);
		
		//update PID in case another process took over
		if (file_exists($pid_file))
		{
			$pid_old = trim(file_get_contents($pid_file));
		}
	}
}
if (!file_put_contents($pid_file, getmypid(), LOCK_EX))
{
	trigger_error("Could not create PID file {$pid_file}", E_USER_ERROR);
}
$pid_mod = substr(sprintf('%o', fileperms($pid_file)), -4);
if ($pid_mod!="0777")
{
	chmod($pid_file, 0777);
}

//copy reference genome and index files
print "Copying genome files...\n";
list($files) = exec2("ls {$genome_folder}/{$build}.*");
foreach($files as $file)
{
	$base = basename($file);
	print "  rsync-ing genome file '$base'.\n";
	$existed_before = file_exists($local_data."/".$base);
	list($stdout, $stderr) = exec2("{$rsync} {$genome_folder}{$base} {$local_data}/{$base}");
	foreach(array_merge($stdout, $stderr) as $line)
	{
		$line = trim($line);
		if ($line=="") continue;
		print "    $line\n";
	}
	if (!$existed_before)
	{
		if (!chmod($local_data."/".$base, 0777))
		{
			trigger_error("Could not change privileges of local data folder '{$local_data}/{$base}'!", E_USER_ERROR);
		}
	}
}

######################### VEP cache #########################
if ($build=="GRCh38")
{
	$annotation_folder = get_path("vep_data")."/cache/";
	$local_annotation_folder = "{$local_data}/".basename(get_path("vep_data"))."/";

	print "\n";
	print "### VEP cache ###\n";
	print "from: {$annotation_folder}\n";
	print "to  : {$local_annotation_folder}\n";
	print "\n";

	//create local annotation folder if missing
	if (!file_exists($local_annotation_folder))
	{
		if (!mkdir($local_annotation_folder))
		{
			trigger_error("Could not create local data annotation folder '{$local_annotation_folder}'!", E_USER_ERROR);
		}
		if (!chmod($local_annotation_folder, 0777))
		{
			trigger_error("Could not change privileges of local annotation data folder '{$local_annotation_folder}'!", E_USER_ERROR);
		}
	}
	
	//remove outdated annotation data
	$update = false;
	$info = "/homo_sapiens/112_GRCh38/info.txt"; //ensembl-vep-112
	if (file_exists("{$local_annotation_folder}/{$info}"))
	{
		exec("diff {$local_annotation_folder}/{$info} {$annotation_folder}/{$info}", $output, $code);
		if ($code==0)
		{
			print "Annotation infos in '$info' match.\n";
		}
		else
		{
			print "Annotation infos in '$info' differ. Deleting old data and performing update!\n";
			exec2("rm -rf {$local_annotation_folder}/*");
			$update = true;
		}
	}
	else
	{
		// info.txt missing -> perform update
		print "{$local_annotation_folder}/{$info} missing. Performing update!\n";
		$update = true;
	}
	
	//check that at least one rsync finished
	$finished = "/rsync_done.txt";
	if (!file_exists("{$local_annotation_folder}/{$finished}"))
	{
		print "No rsync finished. Performing update!\n";
		$update = true;
	}
		
	//perform update
	if ($update)
	{
		print "rsync-ing annotation data...\n";
		
		list($stdout, $stderr) = exec2("{$rsync} {$annotation_folder} {$local_annotation_folder}");
		foreach(array_merge($stdout, $stderr) as $line)
		{
			$line = trim($line);
			if ($line=="") continue;
			
			print "  {$line}\n";
		}
		
		touch("{$local_annotation_folder}/{$finished}");
	}
	else
	{
		print "No rsync necessary!\n";
	}
}

######################### VEP annotation databases #########################
if ($build=="GRCh38")
{
	if (get_path("copy_dbs_to_local_data"))
	{
		$local_annotation_folder = "{$local_data}/ensembl-vep-dbs/";

		//create local folder if missing
		if (!file_exists($local_annotation_folder))
		{
			if (!mkdir($local_annotation_folder))
			{
				trigger_error("Could not create local data annotation folder '{$local_annotation_folder}'!", E_USER_ERROR);
			}
			if (!chmod($local_annotation_folder, 0777))
			{
				trigger_error("Could not change privileges of local annotation data folder '{$local_annotation_folder}'!", E_USER_ERROR);
			}
		}
		
		print "\n";
		print "### VEP annotation databases ###\n";
		print "to  : {$local_annotation_folder}\n";
		print "\n";
			
		print "rsync-ing annotation databases...\n";
		
		//determine databases to sync
		foreach($db_files as $db_file)
		{
			$source = $data_folder.$db_file;
			$target = $local_annotation_folder.basename($db_file);
			print "  rsync-ing database {$source}\n";
			list($stdout, $stderr) = exec2("{$rsync} {$source} {$target}");
			foreach(array_merge($stdout, $stderr) as $line)
			{
				$line = trim($line);
				if ($line=="") continue;
				
				print "    {$line}\n";
			}
			
			$index_format = ".tbi";
			if (file_exists($source.".csi"))
			{
				$index_format = ".csi";
			}

			$source_index = $source.$index_format;

			if (file_exists($source_index))
			{
				print "  rsync-ing database index {$source_index}\n";
				
				$target_index = $target.$index_format;

				list($stdout, $stderr) = exec2("{$rsync} {$source_index} {$target_index}");
				foreach(array_merge($stdout, $stderr) as $line)
				{
					$line = trim($line);
					if ($line=="") continue;
					
					print "    {$line}\n";
				}
			}
		}
	}
}


######################### check DB files #########################
if ($build=="GRCh38" && $check)
{
	print "\n";
	print "### checking annotation databases ###\n";
	$genome = genome_fasta($build);
	
	//determine databases to sync
	foreach($db_files as $db_file)
	{
		$filename = realpath($data_folder.$db_file);
		$basename = basename($filename);
		if (ends_with($filename, ".vcf.gz"))
		{
			$log = "./".$basename."_check.txt";
			print "Checking: $basename ($log)\n";
			$parser->execApptainer("ngs-bits", "VcfCheck", "-in $filename -lines 5000000 -ref $genome > $log 2>&1", [$filename, $genome], [dirname($log)], false, true, false);
			
			$errors = 0;
			$warnings = 0;
			$h = fopen($log, 'r');
			while(!feof($h))
			{
				$line = fgets($h);
				if (starts_with($line, "ERROR:")) ++$errors;
				if (starts_with($line, "WARNING:")) ++$warnings;
			}
			if ($errors>0) print "  Errors: $errors\n";
			if ($warnings>0) print "  Warnings: $warnings\n";
		}
		else
		{
			print "Skipped: $basename\n";
		}
	}
}

######################### copy apptainer containers #########################
if (get_path("copy_dbs_to_local_data"))
{
	$network_folder = get_path("container_folder");
	$network_checksum_folder = $network_folder."/checksums/";
	$local_folder = get_path("local_data")."/container/";
	$local_checksum_folder = $local_folder."/checksums/";

	print "\n";
	print "### Copy apptainer containers ###\n";
	print "from: {$network_folder}\n";
	print "to  : {$local_folder}\n";

	// Check if the network container folder exists
	if (!file_exists($network_folder))
	{
		trigger_error("Container folder not found: {$network_folder}. The Apptainer containers may not have been downloaded yet.", E_USER_ERROR);
	}
	
	if (!file_exists($network_checksum_folder)) //Check if checksum folder exists
	{
		if (!mkdir($network_checksum_folder))
		{
			trigger_error("Could not create checksum folder '{$network_checksum_folder}'!", E_USER_ERROR);
		}
		if (!chmod($network_checksum_folder, 0777))
		{
			trigger_error("Could not change privileges of checksum folder '{$network_checksum_folder}'!", E_USER_ERROR);
		}
	}

	// Create local container folder
	if (!file_exists($local_folder))
	{
		if (!mkdir($local_folder))
		{
			trigger_error("Could not create local data folder '{$local_folder}'!", E_USER_ERROR);
		}
		if (!chmod($local_folder, 0777))
		{
			trigger_error("Could not change privileges of local data folder '{$local_folder}'!", E_USER_ERROR);
		}
	}

	if (!file_exists($local_checksum_folder)) //Check if checksum folder exists
	{
		if (!mkdir($local_checksum_folder))
		{
			trigger_error("Could not create checksum folder '{$local_checksum_folder}'!", E_USER_ERROR);
		}
		if (!chmod($local_checksum_folder, 0777))
		{
			trigger_error("Could not change privileges of checksum folder '{$local_checksum_folder}'!", E_USER_ERROR);
		}
	}

	// get list of apptainer containers to transfer from settings file
	$ini = get_ini();
	foreach($ini as $key => $value)
	{
		if (!starts_with($key, "container_")) continue;
		if ($key=="container_folder") continue;
		
		$container_file = $network_folder."/".substr($key, 10)."_".$value.".sif";
		$base = basename($container_file);
		$local_container_file = $local_folder.$base;
		$md5_local = "{$local_checksum_folder}{$base}.md5";
		$md5_net = "{$network_checksum_folder}{$base}.md5";

		//Check that container exists in network folder
		if (!file_exists($container_file)) trigger_error("Could not find container $container_file. Make sure your settings.ini is correct and the container exists in $network_folder");

		//Make sure that md5 checksum for container in network directory is present
		if (!file_exists($md5_net))
		{
			print "MD5 checksum missing for $container_file, creating it now.\n";
			exec2("md5sum -b {$container_file} | cut -d ' ' -f1 > $md5_net");
		}

		//Make sure md5 checksum for container in local data is present
		if (file_exists($local_container_file) && !file_exists($md5_local))
		{
			print "MD5 checksum missing for $local_container_file, creating it now.\n";
			exec2("md5sum -b {$local_container_file} | cut -d ' ' -f1 > $md5_local");
		}

		//Check md5 sum and update local container folder
		if (file_exists($md5_local))
		{
			exec("diff $md5_local $md5_net", $output, $code);
			if ($code==0)
			{
				print "MD5 checksums of container '{$base}.md5' match - local copy is up-to-date\n";
				continue;
			}
			else
			{
				print "MD5 checksums of container '{$base}.md5' differ. Updating container {$base}!\n";
				// Copy the new container version
				list($stdout, $stderr) = exec2("cp {$container_file} {$local_container_file}");
				foreach (array_merge($stdout, $stderr) as $line) 
				{
					$line = trim($line);
					if ($line == "") continue;
					print "    $line\n";
				}

				// Set permissions on the new local copy
				@chmod($local_container_file, 0777);
				
				//create new local md5 sum
				exec2("rm $md5_local");
				exec2("md5sum -b {$local_container_file} | cut -d ' ' -f1 > $md5_local");
			}
		}
		else
		{
			print "Container $base missing in local data folder. Copying container {$base}!\n";
			// Copy the new container
			list($stdout, $stderr) = exec2("cp {$container_file} {$local_container_file}");
			foreach (array_merge($stdout, $stderr) as $line) 
			{
				$line = trim($line);
				if ($line == "") continue;
				print "    $line\n";
			}

			// Set permissions on the new local copy
			@chmod($local_container_file, 0777);
			
			//create local md5 sum
			exec2("md5sum -b {$local_container_file} | cut -d ' ' -f1 > $md5_local");
		}
	}
}

######################### Remove PID file if it is still ours #########################
$current_pidfile_pid = trim(file_get_contents($pid_file));
if ($current_pidfile_pid==getmypid())
{
	unlink($pid_file);
}
else
{
	print "PID file was overwritten by another process ($current_pidfile_pid). Not removing it.\n";
}

?>
