<?php
/** 
	@page index_genome
*/

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

//parse command line arguments
$parser = new ToolBase("index_genome", "Indexes a genome FASTA file with BWA and samtools.");
$parser->addInfile("in", "FASTA file to index.", false);
$parser->addFlag("mask", "Mask false duplications in GRCh38 (see https://www.nature.com/articles/s41587-021-01158-1).");
extract($parser->parse($argv));

if ($mask)
{
	$exclusion_bed = $parser->tempFile("_exclusion.bed");
	exec2("wget -O {$exclusion_bed} https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_GRC_exclusions.bed");
	$tmp = $parser->tempFile("_masked.fa");
	$parser->execApptainer("ngs-bits", "FastaMask", "-in {$in} -reg {$exclusion_bed} -out {$tmp}", [$in]);
	$parser->moveFile($tmp, $in);
}

//BWA index
if (get_path("use_bwa1"))
{
	$parser->execApptainer("bwa", "bwa index", "-a bwtsw {$in}", [$in]);
}
else
{
	$suffix = trim(get_path("bwa_mem2_suffix", false));
	$parser->execApptainer("bwa-mem2", "bwa-mem2{$suffix}", " index {$in}", [$in]);
}

//samtools FAI file
$parser->execApptainer("samtools", "samtools", "faidx {$in}", [$in]);
exec2("md5sum -b {$in} > {$in}.md5");

//GATK dict file
$parser->execApptainer("gatk", "gatk", "CreateSequenceDictionary -R {$in}", [$in]);

//create samtools ref_cache
$parser->execApptainer("samtools", "seq_cache_populate.pl", "-root ".dirname($in)."/samtools_ref_cache {$in}", [$in]);

?>
