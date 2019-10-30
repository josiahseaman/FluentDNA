<?php
ini_set('max_execution_time', '300');

//Open sequence.fasta and read it
//get start and end as parameters
//ignore the > heading
//count GC percentage from start to end
//return percentage

$start=0;
$end=0;
$filepath="";
$error=0;

if(isset($_GET["start"])) {$start=$_GET["start"];}
if(isset($_GET["end"])) {$end=$_GET["end"];}
if(isset($_GET["filepath"])) {$filepath=$_GET["filepath"];}

$filepath = preg_replace ( "/.*(?=dnadata)/" , "" , $filepath );

$file=$filepath."sequence.fasta";

$g=0;
$c=0;
$t=0;
$a=0;
$n=0;
$total=0;

if ($filepath) {$handle = fopen($file, "r");}
if ($handle) {
    while (($line = fgets($handle)) !== false) {
    	if (substr($line, 0, 1) == '>') { continue; }
    	if ($total >= $end) { break; /*done*/}
        // process the line read.
        $line=trim($line);
        
        			//add the values in the range only in the row 
        			//for all the characters in the row
        			for ($i = 0, $j = strlen($line); $i < $j; $i++) {
        				  if ($total < $start) {$total++;continue;}
        				  
     							$a = $a + substr_count($line[$i], 'A'); 
			        		$c = $c + substr_count($line[$i], 'C'); 
			        		$g = $g + substr_count($line[$i], 'G'); 
			        		$t = $t + substr_count($line[$i], 'T'); 
			        		$n = $n + substr_count($line[$i], 'N'); 
			        		
			        		if ($total >= $end) {break;}
			        		$total++;
							}				
			       
      	
    }
} else {
    echo("Error: couldn't open sequence");
    echo('<div>Sequence: '.$file.'</div>');
    $error=1;
} 

fclose($handle);

$allCountedNucleotides = $a+$g+$c+$t+$n;

$gc = ($g+$c)/$allCountedNucleotides;
$at = ($a+$t)/$allCountedNucleotides;

$pa = $a/$allCountedNucleotides;
$pc = $c/$allCountedNucleotides;
$pt = $t/$allCountedNucleotides;
$pg = $g/$allCountedNucleotides;
$pn = $n/$allCountedNucleotides;

if (!$error){
echo('<div>Nucleic density (start = '.$start.' , end = '.$end.')</div>');
echo('<table style="width:600px;"><tr><td>A</td><td>C</td><td>T</td><td>G</td><td>N</td><td>G+C</td><td>A+T</td></tr>');
echo("<tr><td>".round($pa,3)."</td><td>".round($pc,3)."</td><td>".round($pt,3)."</td><td>".round($pg,3)."</td><td>".round($pn,3)."</td><td>".round($gc,3)."</td><td>".round($at,3));
echo('</table>');
}

?>