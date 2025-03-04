<?php

require_once("framework.php");

$name = "an_scarHRD";
start_test($name);

//init NGSD
init_ngsd($name);

//empty
$log_file = output_folder()."/{$name}_out0.log";
check_exec("php ".src_folder()."/Tools/{$name}.php -cnvs ".data_folder()."/{$name}_in0.tsv -tumor DNA220000_01 -normal DNA220000_01 -out_folder ".output_folder()."/ -db NGSD_TEST --log $log_file");
check_file(output_folder()."/DNA220000_01-DNA220000_01_HRDresults.txt", data_folder().$name."_out0.txt");

// all cnvs
$log_file = output_folder()."/{$name}_out1.log";
check_exec("php ".src_folder()."/Tools/{$name}.php -skip_prefilter -cnvs ".data_folder()."/{$name}_in1.tsv -tumor DNA220001_01 -normal DNA220002_01 -out_folder ".output_folder()."/ -db NGSD_TEST --log $log_file");
check_file(output_folder()."/DNA220001_01-DNA220002_01_HRDresults.txt", data_folder().$name."_out1.txt");

// report filtered - remove cnvs that are marked as artifacts in NGSD
$log_file = output_folder()."/{$name}_out2.log";
check_exec("php ".src_folder()."/Tools/{$name}.php -cnvs ".data_folder()."/{$name}_in2.tsv -tumor DNA220003_01 -normal DNA220004_01 -out_folder ".output_folder()."/ -db NGSD_TEST -report_filtered --log $log_file");
check_file(output_folder()."/DNA220003_01-DNA220004_01_HRDresults.txt", data_folder().$name."_out2.txt");

// prefilter cnvs
$log_file = output_folder()."/{$name}_out3.log";
check_exec("php ".src_folder()."/Tools/{$name}.php -cnvs ".data_folder()."/{$name}_in3.tsv -tumor DNA220011_01 -normal DNA220022_01 -out_folder ".output_folder()."/ -db NGSD_TEST --log $log_file");
check_file(output_folder()."/DNA220011_01-DNA220022_01_HRDresults.txt", data_folder().$name."_out3.txt");

end_test();

?>