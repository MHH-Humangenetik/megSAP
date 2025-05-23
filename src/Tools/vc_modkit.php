<?php 
/** 
	@page vc_modkit
*/

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

// parse command line arguments
$parser = new ToolBase("vc_modkit", "ONT methylation annotation using modkit.");
$parser->addInfile("bam",  "Input file in BAM format. Note: .bam.bai file is required!", false);
$parser->addOutfile("bed", "Output BED file containing the methylation info. Bgzipped and indexed.", false);

//optional
$parser->addOutfile("summary", "Optional summary file.", true);
$parser->addInt("threads", "The maximum number of threads used.", true, 1);
$parser->addString("build", "The genome build to use.", true, "GRCh38");

extract($parser->parse($argv));

//init
$genome = genome_fasta($build);
$log_file = $parser->tempFile("_modkit_pileup.log");
$uncompressed_bed = $parser->tempFile("_modkit.bed");

//run modkit
$args = [];
$args[] = "pileup";
$args[] = $bam;
$args[] = $uncompressed_bed;
$args[] = "--ref ".$genome;
$args[] = "--log-filepath ".$log_file;
$args[] = "--threads ".$threads;

//set bind paths for modkit
$in_files = array();
$in_files[] = $bam;
$in_files[] = $genome;

$parser->execApptainer("modkit", "modkit", implode(" ", $args), $in_files);

//copy logfile 
$parser->log("modkit pileup log file", file($log_file));

//store compressed file
$parser->execApptainer("htslib", "bgzip", "-c $uncompressed_bed > $bed", [], [dirname($bed)]);
$parser->execApptainer("htslib", "tabix", "-f -p bed $bed", [], [dirname($bed)]);

//run summary
if (isset($summary))
{
	$log_file2 = $parser->tempFile("_modkit_summary.log");
	$args = [];
	$args[] = "summary";
	$args[] = $bam;
	$args[] = "--threads ".$threads;
	$args[] = "--log-filepath ".$log_file2;
	$args[] = " > ".$summary;

	$parser->execApptainer("modkit", "modkit", implode(" ", $args), $in_files, [dirname($summary)]);
	//copy logfile 
	$parser->log("modkit summary log file", file($log_file2));
}



?>
