#!/usr/bin/env perl

# NTM
# 27/07/2020


# This is a wrapper script for the grexome-TIMC secondary analysis
# pipeline.
# Args: see $USAGE.


use strict;
use warnings;
use Getopt::Long;
use POSIX qw(strftime);
use File::Copy qw(copy);
use File::Basename qw(basename);
use File::Temp qw(tempdir);
use FindBin qw($RealBin);

use lib "$RealBin";
use grexome_metaParse qw(parsePathologies parseSamples parseCandidateGenes);


# we use $0 in every stderr message but we really only want
# the program name, not the path
$0 = basename($0);


#############################################
## hard-coded stuff that shouldn't change much

# number of parallel jobs to run for the initial "bgzip -d",
# and then for steps steps 1 and 6.
# These defaults are good for us (dual Xeon 4114) but tweaking
# could improve performance depending on your hardware.
my $numJobsGunzip = 6;
my $numJobs1 = 20;
my $numJobs6 = 16;

# name(+path if needed) of gzip-like binary, we like bgzip (multi-threaded)
my $bgzip = "bgzip";
(`which $bgzip` =~ /$bgzip/) ||
    die "E $0: the bgzip executable $bgzip can't be found\n";
# full command. If you don't have bgzip you could use gzip but without --threads
$bgzip = "$bgzip -cd --threads $numJobsGunzip ";


#############################################
## options / params from the command-line

# metadata file with all samples
my $samples;

# pathologies metadata file
my $pathologies;

# comma-separated list of path+files holding known candidate genes
my $candidateGenes;

# input bgzipped multi-sample GVCF or VCF
my $inFile;

# input multi-sample VCF containing CNV calls, possibly (b)gzipped
my $cnvs;

# input Interactome Results File(.tsv) file produced by InteractomeResults.py
my $ininteractomeresults;

# outDir must not exist, it will be created and populated with
# subdirs (containing the pipeline results), logfiles (in debug mode),
# and copies of all provided metadata files.
my $outDir;

# path+file of the config file holding all install-specific params,
# defaults to the distribution-povided file that you can edit but
# you can also copy it elsewhere and customize it, then use --config
my $config = "$RealBin/grexomeTIMCsec_config.pm";

# if $canon is true, results are restricted to CANONICAL transcripts.
# Otherwise you can filter the CSVs using the CANONICAL column,
# and we also still produce canonical-only Cohorts files (because
# the unrestricted Cohorts files can be too heavy for oocalc/excel)
my $canon = '';

# debug: if true:
# - run each step one after the other (no pipes)
# - check result of previous step before starting next (die if prev step failed)
# - keep all intermediate files (no cleanup)
# - produce individual logfiles for each step
my $debug = '';

# $debugVep: activate --debug specifically for 3_runVep.pl
my $debugVep = '';

# help: if true just print $USAGE and exit
my $help = '';

my $USAGE = "Parse a GVCF or VCF, run the complete grexome-TIMC secondary analysis pipeline,
and produce results (+ logs and copies of the metadata) in the provided outDir (which must not exist).
Each step of the pipeline is a stand-alone self-documented script, this is just a wrapper.
Every install-specific param (eg paths to shared data) should be in grexomeTIMCsec_config.pm.
Arguments [defaults] (all can be abbreviated to shortest unambiguous prefixes):
--samples : samples metadata xlsx file, with path
--pathologies : [optional] pathologies metadata xlsx file, with path
--candidateGenes : [optional] known candidate genes in xlsx files, comma-separated, with paths
--infile : bgzipped multi-sample GVCF or VCF file to parse
--cnvs : [optional] multi-sample VCF file containing CNV calls, possibly (b)gzipped, conventions are:
       	 ALT is <DEL> or <DUP>, INFO contains END=, FORMAT must start with GT:RR and contain BF
--ininteractomeresults : [optional] Interactome Results File(.tsv) file produced by InteractomeResults.py - [For integrating Interactome data]
--outdir : subdir where results will be created, must not pre-exist
--config [defaults to grexomeTIMCsec_config.pm alongside this script] : your customized copy (with path) of the distributed *config.pm
--canonical : restrict results to canonical transcripts
--debug : activate debug mode => slower, keeps all intermediate files, produce individual logfiles
--debugVep : activate debug mode specifically for VEP => compare VEP annotations between runs
--help : print this USAGE";

GetOptions ("samples=s" => \$samples,
            "pathologies=s" => \$pathologies,
            "candidateGenes=s" => \$candidateGenes,
            "infile=s" => \$inFile,
            "cnvs=s" => \$cnvs,
            "ininteractomeresults=s" => \$ininteractomeresults,
            "outdir=s" => \$outDir,
            "config=s" => \$config,
            "canonical" => \$canon,
            "debug" => \$debug,
            "debugVep" => \$debugVep,
            "help" => \$help)
    or die("E $0: Error in command line arguments\n$USAGE\n");

# make sure required options were provided and sanity check them
($help) && die "$USAGE\n\n";

($samples) || die "E $0: you must provide a samples file\n\n$USAGE\n";
(-f $samples) || die "E $0: the supplied samples file doesn't exist:\n$samples\n";

($inFile) || die "E $0: you must provide an input bgzipped (G)VCF file\n";
(-f $inFile) || die "E $0: the supplied infile doesn't exist\n";
($inFile =~ /\.gz$/) || die "E $0: the supplied infile doesn't seem bgzipped\n";

if ($cnvs) {
    (-f $cnvs) || die "E $0: the supplied cnvs file doesn't exist\n";
}

if ($ininteractomeresults) {
    (-f $ininteractomeresults) || die "E $0: the supplied interactome results file doesn't exist\n";
}

# immediately import $config, so we die if file is broken
# if $config doesn't have a path component, prepend ./ to avoid loading the dist version
# (in case the dist was copied into current dir and customized but not renamed)
($config =~ m~/~) || ($config = "./$config");
(-f $config) ||  die "E $0: the supplied config.pm doesn't exist: $config\n";
require($config);
grexomeTIMCsec_config->import(qw(refGenome vepCacheFile vepPluginDataPath fastTmpPath),
			      qw(coveragePath gtexDatafile gtexFavoriteTissues subCohorts));

($outDir) || die "E $0: you must provide an outDir\n";
(-e $outDir) &&
    die "E $0: outDir $outDir already exists, remove it or choose another name.\n";
mkdir($outDir) || die "E $0: cannot mkdir outDir $outDir\n";


# copy all provided metadata files into $outDir
copy($samples, $outDir) ||
    die "E $0: cannot copy samples metadata to outDir: $!\n";
# use the copied versions in scripts (eg if original gets edited while analysis is running...)
$samples = "$outDir/".basename($samples);

if ($pathologies) {
    (-f $pathologies) || die "E $0: the supplied pathologies file doesn't exist\n";
    copy($pathologies, $outDir) ||
	    die "E $0: cannot copy pathologies metadata to outDir: $!\n";
    $pathologies = "$outDir/".basename($pathologies);
}

my @candNew = ();
if ($candidateGenes) {
    foreach my $candFile (split(/,/, $candidateGenes)) {
	(-f $candFile) ||
	    die "E $0: the supplied candidateGenes file $candFile doesn't exist\n";
	copy($candFile, $outDir) ||
	    die "E $0: cannot copy candidateGenes file $candFile to outDir: $!\n";
	# use the copies in script
	push(@candNew, "$outDir/".basename($candFile));
    }
    $candidateGenes = join(',', @candNew);
}

my $now = strftime("%F %T", localtime);
warn "I $now: $0 - starting to run\n";

# sanity-check all provided metadata files: just parse and ignore the results,
# we'll die if anything is wrong 
if ($pathologies) {
    &parsePathologies($pathologies);
    &parseSamples($samples, $pathologies);
    if ($candidateGenes) {
	&parseCandidateGenes($candidateGenes, $samples, $pathologies);
    }
}
else {
    &parseSamples($samples);
    if ($candidateGenes) {
	&parseCandidateGenes($candidateGenes, $samples);
    }
}

# sanity-check: any sample listed in $samplesFile MUST be present in $inFile,
# otherwise we can get errors in later steps (eg step 7).
# We also need the number of samples in $samples, to set $min_hr (for filtering).
my $numSamples = 0;
{
    # grab the list of samples present in $inFile:
    my @samplesFromInFile = split(/\s+/, `zgrep --max-count=1 '#CHROM' $inFile`);
    # first 9 columns aren't sampleIDs
    splice(@samplesFromInFile,0,9);
    # key == sampleID present in $inFile, value==1
    my %samplesFromIn = ();
    foreach my $sample (@samplesFromInFile) {
	$samplesFromIn{$sample} = 1;
    }

    # grab the samples listed in $samplesFile: just use the first hashref from
    # parseSamples, keys are samples
    my ($samplesFromMetadataR) = &parseSamples($samples);
    $numSamples = scalar(keys(%$samplesFromMetadataR));
    foreach my $s (sort keys(%$samplesFromMetadataR)) {
	(defined $samplesFromIn{$s}) ||
	    die "E $0: every sample defined in samples file $samples MUST be present in inFile $inFile, sample $s isn't.\n";
    }
}

# variant-caller string/name, will be added to all final filenames
my $caller;
# try to find the caller name in $inFile
foreach my $c ("Strelka", "GATK", "ElPrep", "DeepVariant") {
    if ($inFile =~ /$c/i) {
	$caller = $c;
	warn "I $0: variant-caller id $caller will be appended to all filenames\n";
	last;
    }
}
($caller) ||
    warn "W $0: variant-caller could not be auto-detected, filenames won't be tagged\n";


#############################################

# all intermediate results / tmp working files are created in $tmpdir,
# a randomly-named subdir of &fastTmpPath() to avoid clashes
# $tmpdir is removed afterwards except in debug mode
my $tmpdir = tempdir(DIR => &fastTmpPath());

######################
# STEPS 1-6, piped into each other except in debug mode

# decompress infile and step 1
my $com = "$bgzip $inFile | perl $RealBin/1_filterBadCalls.pl --samplesFile=$samples --tmpdir=$tmpdir/FilterTmp/ --jobs $numJobs1 ";
if ($debug) {
    # specific logfile from step and save its output
    $com .= "2> $outDir/step1.err > $outDir/step1.out";
    system($com) && die "E $0: debug mode on, step1 failed: $?";
    # next step will read this step's output
    $com = "cat $outDir/step1.out ";
}

# step 1B if we have a CNVs file
if ($cnvs) {
    $com .= " | perl $RealBin/1B_collateVCFs.pl --vcf=$cnvs ";
    if ($debug) {
	$com .= "2> $outDir/step1B.err > $outDir/step1B.out";
	system($com) && die "E $0: debug mode on, step1B failed: $?";
	$com = "cat $outDir/step1B.out ";
    }
}

# step 2
$com .= " | perl $RealBin/2_sampleData2genotypes.pl ";
if ($debug) {
    $com .= "2> $outDir/step2.err > $outDir/step2.out";
    system($com) && die "E $0: debug mode on, step2 failed: $?";
    $com = "cat $outDir/step2.out ";
}

# step 3
$com .= " | perl $RealBin/3_runVEP.pl --cacheFile=".&vepCacheFile()." --genome=".&refGenome()." --dataDir=".&vepPluginDataPath()." --tmpDir=$tmpdir/runVepTmpDir/ ";
($debugVep) && ($com .= "--debug ");
if ($debug) {
    $com .= "2> $outDir/step3.err > $outDir/step3.out";
    system($com) && die "E $0: debug mode on, step3 failed: $?";
    $com = "cat $outDir/step3.out ";
}

# step 4
$com .= " | perl $RealBin/4_vcf2tsv.pl ";
if ($debug) {
    $com .= "2> $outDir/step4.err > $outDir/step4.out";
    system($com) && die "E $0: debug mode on, step4 failed: $?";
    $com = "cat $outDir/step4.out ";
}

# step 5.1
$com .= " | perl $RealBin/5.1_addGTEX.pl --favoriteTissues=".&gtexFavoriteTissues()." --gtex=".&gtexDatafile($RealBin)." ";
if ($debug) {
    $com .= "2> $outDir/step5.1.err > $outDir/step5.1.out";
    system($com) && die "E $0: debug mode on, step5.1 failed: $?";
    $com = "cat $outDir/step5.1.out ";
}

# step 5.2
$com .= " | python3 $RealBin/5.2_addInteractome.py --ininteractomeresults $ininteractomeresults ";
if ($debug) {
    $com .= "2> $outDir/step5.2.err > $outDir/step5.2.out";
    system($com) && die "E $0: debug mode on, step5.2 failed: $?";
    $com = "cat $outDir/step5.2.out ";
}

# step 6
$com .= " | perl $RealBin/6_extractCohorts.pl --samples=$samples ";
($pathologies) && ($com .= "--pathologies=$pathologies ");
($candidateGenes) && ($com .= "--candidateGenes=$candidateGenes ");
$com .= "--outDir=$tmpdir/Cohorts/ --tmpDir=$tmpdir/TmpExtract/ --jobs=$numJobs6 ";
if ($debug) {
    $com .= "2> $outDir/step6.err";
    system($com) && die "E $0: debug mode on, step6 failed: $?";
}
else {
    # after step6 we have several files (one per cohort) so no more piping,
    # run steps 1-6, all logs go to stderr
    system($com) && die "E $0: steps 1 to 6 seem to have failed? $?";
}

######################
# subsequent steps work on the individual CohortFiles

# STEP 7: filter variants and reorder columns, clean up unfiltered CohortFiles
$com = "perl $RealBin/7_filterAndReorderAll.pl --indir $tmpdir/Cohorts/ --outdir $tmpdir/Cohorts_Filtered/ --reorder ";
# set min_hr to 20% of $numSamples
my $min_hr = int($numSamples * 0.2);
$com .= "--min_hr=$min_hr ";
# other hard-coded filters:
$com .= "--max_ctrl_hv 3 --max_ctrl_het 10 --no_mod ";
$com .= " --max_af_gnomad 0.01 --max_af_1kg 0.03 ";
($canon) && ($com .= "--canonical ");
if ($debug) {
    $com .= "2> $outDir/step7.err";
}
else {
    # remove unfiltered results in non-debug mode
    $com .= " ; rm -r $tmpdir/Cohorts/";
}
system($com) && die "E $0: step7 failed: $?";

# STEP 8 - SAMPLES
$com = "perl $RealBin/8_extractSamples.pl --samples $samples ";
$com .= "--indir $tmpdir/Cohorts_Filtered/ --outdir $outDir/Samples/ ";
(&coveragePath()) && ($com .= "--covdir ".&coveragePath()." ");
if ($debug) {
    $com .= "2> $outDir/step8s.err";
}
system($com) && die "E $0: step8-samples failed: $?";

# STEP 8 - TRANSCRIPTS , adding patientIDs
$com = "( perl $RealBin/8_extractTranscripts.pl --indir $tmpdir/Cohorts_Filtered/ --outdir $tmpdir/Transcripts_noIDs/ ";
($pathologies) && ($com .= "--pathologies=$pathologies ");
$com .= "; perl $RealBin/8_addPatientIDs.pl $samples $tmpdir/Transcripts_noIDs/ $outDir/Transcripts/ )";
if ($debug) {
    $com .= "2> $outDir/step8t.err";
}
else {
    $com .= " ; rm -r $tmpdir/Transcripts_noIDs/ ";
}
system($com) && die "E $0: step8-transcripts failed: $?";

# STEP 9 - FINAL COHORTFILES , require at least one HV or HET sample and add patientIDs
$com = "( perl $RealBin/9_requireUndiagnosed.pl $tmpdir/Cohorts_Filtered/ $tmpdir/Cohorts_FINAL/ ; ";
$com .= " perl $RealBin/8_addPatientIDs.pl $samples $tmpdir/Cohorts_FINAL/ $outDir/Cohorts/ )";
if ($debug) {
    $com .= "2> $outDir/step9-finalCohorts.err";
}
else {
    $com .= " ; rm -r $tmpdir/Cohorts_Filtered/ $tmpdir/Cohorts_FINAL/ ";
}
system($com) && die "E $0: step9-finalCohorts failed: $?";

# STEP 9 - COHORTFILES_CANONICAL , if called without --canon we still
# produce CohortFiles restricted to canonical transcripts (in case
# the AllTranscripts CohortFiles are too heavy for oocalc/excel)
if (! $canon) {
    $com = "perl $RealBin/7_filterAndReorderAll.pl --indir $outDir/Cohorts/ --outdir $outDir/Cohorts_Canonical/ --canon ";
    if ($debug) {
	$com .= "2> $outDir/step9-cohortsCanonical.err";
    }
    system($com) && die "E $0: step9-cohortsCanonical failed: $?";
}

######################
# STEP 9 - SUBCOHORTS (can run after requireUndiagnosed and addPatientIDs)
#
# Only runs if sub-cohorts are defined in &subCohorts() (and the files exist):
# the idea is to produce Cohorts and Transcripts files corresponding to subsets of
# samples, all affected with the same pathology. Typically the subsets are samples
# that were provided by collaborators, and this allows us to send them the results
# concerning their patients.

# key==path+file defining a subCohort, value==pathologyID
my $subCohortsR = &subCohorts();

# don't do anything if no subcohort file exists
my $doSubCs = 0;
foreach my $subC (keys(%$subCohortsR)) {
    if (-e $subC) {
	$doSubCs=1;
    }
    else {
	warn "W $0: sub-cohort file $subC defined in \&subCohorts() but this file doesn't exist. Skipping this sub-cohort.\n";
    }
}
if ($doSubCs) {
    mkdir("$outDir/SubCohorts") || die "E $0: cannot mkdir $outDir/SubCohorts\n";
    $com = "( ";
    foreach my $subC (keys(%$subCohortsR)) {
	my $patho = $subCohortsR->{$subC};
	# grab filename from $subC and remove leading "subCohort_" and trailing .txt
	my $outFileRoot = basename($subC);
	($outFileRoot =~ s/^subCohort_//) || die "E $0: cannot remove leading subCohort_ from subCohortFile $outFileRoot\n";
	($outFileRoot =~ s/\.txt$//) || die "E $0: cannot remove .txt from subCohortFile $outFileRoot\n";
	$outFileRoot = "$outDir/SubCohorts/$outFileRoot";

	$com .= "perl $RealBin/9_extractSubcohort.pl $subC < $outDir/Cohorts/$patho.final.patientIDs.csv > $outFileRoot.cohort.csv ; ";
	($canon) ||
	    ($com .= "perl $RealBin/9_extractSubcohort.pl $subC < $outDir/Cohorts_Canonical/$patho.final.patientIDs.canon.csv > $outFileRoot.cohort.canon.csv ; ");
	$com .= "perl $RealBin/9_extractSubcohort.pl $subC < $outDir/Transcripts/$patho.Transcripts.patientIDs.csv > $outFileRoot.transcripts.csv ; ";
    }
    if ($debug) {
	$com .= "2> $outDir/step9-subCohorts.err )";
    }
    else {
	$com .= " )";
    }
    system($com) && die "E $0: step9-subCohorts failed: $?";
}
else {
    warn "I $0: no existing sub-cohort, step9-subCohorts skipped\n";
}

######################
# STEP10 - clean up and QC

# APPENDVC: append $caller (if it was auto-detected) to all final filenames
if ($caller) {
    open (FILES, "find $outDir/ -name \'*csv\' |") ||
	    die "E $0: step10-appendVC cannot find final csv files with find\n";
    while (my $f = <FILES>) {
	chomp($f);
	my $new = $f;
	($new =~ s/\.csv$/.$caller.csv/) ||
	    die "E $0: step10-appendVC cannot add $caller as suffix to $new\n";
	(-e $new) &&
	    die "E $0: step10-appendVC want to rename to new $new but it already exists?!\n";
	rename($f,$new) ||
	    die "E $0: step10-appendVC cannot rename $f $new\n";
    }
    close(FILES);
}


# REMOVEEMPTY: remove files with no data lines (eg if infile concerned only some samples)
open (FILES, "find $outDir/ -name \'*csv\' |") ||
        die "E $0: step10-removeEmpty cannot find final csv files with find\n";
while (my $f = <FILES>) {
    chomp($f);
    my $wc = `cat $f | wc -l`;
    # there's always 1 header line
    ($wc > 1) || unlink($f) ||
	    die "E $0: step10-removeEmpty cannot unlink $f: $!\n";
}
close(FILES);


# QC_CHECKCAUSAL: report coverage of known causal genes by (severe) variants
# QC report will be printed to $qc_causal file
my $qc_causal = "$outDir/qc_causal.csv";

$com = "perl $RealBin/10_qc_checkCausal.pl --samplesFile=$samples --indir=$outDir/Samples/ ";
$com .= "> $qc_causal";
system($com) && die "E $0: step10-qc_causal failed: $?";


######################
# all done, clean up tmpdir
($debug) || rmdir($tmpdir) || die "E $0: all done but can't rmdir(tmpdir), why?\n";

$now = strftime("%F %T", localtime);
warn "I $now: $0 - ALL DONE, completed successfully!\n";