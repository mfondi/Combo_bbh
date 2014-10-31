#!/usr/bin/perl -I /home/marco/Perl/site/lib -I /home/ema/Perl/site/lib
#ARGUMENTS:
#input_file_1: sequenze di cui si vogliono gli ortologhi
#input_file_2: genoma di provenienza delle sequenze in input_file_1
#input_file_3: genomi in cui cercare gli ortologhi

#OUTPUT
#file BBH_of_query_N.txt.txt con best bi-directional hits per ciascuna query
#altri file di output: le sequenze del file di partenza divise in altrettanti file  numerati in base all'ordine in cui appaiono.
#output riassuntivo tabulare. 
use locale;
#use strict;

# no warnings 'redefine';
# no warnings 'uninitialized';
use lib "./marcomod";
#use genomes_melter;

use Bio::SeqIO;
use Bio::Perl;
use Bio::Tools::Run::StandAloneBlast;
use Bio::Seq;
use Bio::Tools::Blast;
use Bio::DB::GenBank;
use Bio::DB::WebDBSeqI;

my @genomes_list;
my $line;
my $line1;
my $line2;
my $line3;
my $line_first_seq;
my @reverse_blast;
my $line_1;
my $putative_orth;
my @fasta_cmd2;
my $n;
my $i;
my $j;
my @first_sequence;
my $seq_number=0;
my @query_list;
my @seq;
my $state;
my $gi;
my $gi2;
my @first_blast;
my @fasta_cmd;
my @gi_first_blast;
my $gi_query;
my $evalue2;
my $all_genomes_to_probe;
my @gi_BBHs;
my $number_of_genomes;
my @genome_list_n;
print "\n\n******************           bbh_3.0           ******************\n";

my $file_input_1= $ARGV[0];
chomp ($file_input_1);
my $file_input_2= $ARGV[1];
chomp ($file_input_2);

my $file_input_3= $ARGV[2];
chomp ($file_input_3);


if (($ARGV[1]eq "") || ($ARGV[0] eq "") || ($ARGV[2] eq "") ){

print "\n\nError: one ore more arguments are missing:\n\nUsage: perl bbh_3.0.pl <seeds file> <reference genome> <path to genomes folder>\n\n";

exit
}

print "\nCleaning up old files (don't mind about eventual bash errors here!) \n";
system("cp $file_input_1 queries.sequences");
system ("mv TABULAR_FBH_OUTPUT.xls TABULAR_FBH_OUTPUT.xls.old");
system ("mv TABULAR_BBH_OUTPUT.xls TABULAR_BFBH_OUTPUT.xls.old");
system ("mv TABULAR_BBH_OUTPUT.csv TABULAR_BBH_OUTPUT.csv.old");
system("rm *.faa");



print "...OK\n";

system("cp $file_input_1 queries.sequences");
system("ls $ARGV[2] > lista_genomi_2.0");
system("cp $ARGV[2]/*.faa .");
my $file_input_3='lista_genomi_2.0';

open(FASTAREADER2, $file_input_3);
@genomes_list = <FASTAREADER2>;



#my $list_finta=genomes_melter::one_genome_file(@genomes_list);

$number_of_genomes=@genomes_list;
open(FASTAREADER2, "queries.sequences");
@query_list = <FASTAREADER2>;
$n=@query_list;


#foreach  $line(@genomes_list){
#chomp($line);
#$line="$line.fasta"
#}


system("formatdb -i \"$file_input_2\" -n bbh_DB_query -p T -e T -o T");

#considera una sequenza alla volta in un file di seq fasta multiplo, tutte in file separati
########sequence n.1##########
for ($i=0; $i<=$n-1; $i++)  {

	if ($query_list[$i] =~ /^>/)    {
	#print"sequence_found\n";
	$query_list[$i] =~/\|(.+?)\|/;
	chomp($query_list[$i]);
	open (SCRIVI2, ">>TABULAR_BBH_OUTPUT.xls");
	print SCRIVI2 "\t$query_list[$i]";
	close SCRIVI2;
	
	chomp($query_list[$i]);
	open (SCRIVI2, ">>TABULAR_BBH_OUTPUT.csv");
	print SCRIVI2 ";$query_list[$i]";
	close SCRIVI2;

	open (SCRIVI2, ">>TABULAR_FBH_OUTPUT.xls");
	print SCRIVI2 "\t$query_list[$i]";
	close SCRIVI2;
		
	

	$seq_number=$seq_number+1;
	$state=0;
	push(@seq, "$query_list[$i]");
	open (SCRIVI, ">query_seq_$seq_number.txt");
	print SCRIVI @query_list[$i];
	close SCRIVI;
	

		for ($j=$i+1; $j<=$n; $j++){
		
			if(($state==0)&($query_list[$j] =~ /^[A-Z]/)){
			
			$line_1=$line_1.$query_list[$j];
			}
			
			else{$state=1;
			open (SCRIVI, ">>query_seq_$seq_number.txt");
			print SCRIVI "$line_1";
			close SCRIVI; 
			
			$line_1='';
			
			}
		}
	}
} 
print "\n\nOk. Data is stored:\n\n";
print "----SEEDS & DATABASE-------------------\n";
print"number of total sequence to probe is $seq_number\n";
print"number of total genomes  to probe is $number_of_genomes\n";
print "---------------------------------------\n";


#open (SCRIVI2, ">>TABULAR_BBH_OUTPUT.xls");
#print SCRIVI2 "\t";
#close SCRIVI2;


foreach $line(@genomes_list){
chomp($line);
#print "$line\n";
system("formatdb -i \"$line\" -n bbh_DB_$line -p T -e T -o T");
}

foreach $line(@genomes_list){
chomp($line);
push(@genome_list_n, "$line\n");
}

#open (SCRIVI2, ">>TABULAR_BBH_OUTPUT.xls");
#print SCRIVI2 "@genome_list_n\t";
#close SCRIVI2;

#open (SCRIVI2, ">>TABULAR_BBH_OUTPUT.xls");
#print SCRIVI2 "\t";
#close SCRIVI2;

foreach $line(@genomes_list){
chomp($line);
$all_genomes_to_probe="$all_genomes_to_probe $line";
}
system("formatdb -i \"$all_genomes_to_probe\" -n all_genomes_to_probe -p T -e T -o T");








open (SCRIVI2, ">>TABULAR_BBH_OUTPUT.xls");
print SCRIVI2 "QUERY $gi_query\n";
close SCRIVI2;

open (SCRIVI2, ">>TABULAR_BBH_OUTPUT.csv");
print SCRIVI2 "$gi_query\n";
close SCRIVI2;

open (SCRIVI2, ">>TABULAR_FBH_OUTPUT.xls");
print SCRIVI2 "QUERY $gi_query\n";
close SCRIVI2;


#print "\n\n\nROUND $i/$seq_number\n>>> >>> >>> >>> >>> >>>   now probing seed: $gi_query   <<< <<< <<< <<< <<< <<< <<< <<< <<< \n";


foreach $line(@genomes_list)  {

print "\n***** probing genome: $line";
print "			*****  ok\n";


open (SCRIVI2, ">>TABULAR_BBH_OUTPUT.xls");
print SCRIVI2 "$line\t";
close SCRIVI2;
open (SCRIVI2, ">>TABULAR_BBH_OUTPUT.csv");
print SCRIVI2 "$line;";
close SCRIVI2;
open (SCRIVI2, ">>TABULAR_FBH_OUTPUT.xls");
print SCRIVI2 "$line\t";
close SCRIVI2;


	for ($i=1; $i<=$seq_number; $i++)  {
	unless (open(BLASTREADER2, "query_seq_$i.txt") ) {
	print "cannot open query_seq_$i.txt\n";
	exit;
	}
	@first_sequence = <BLASTREADER2>;
	close BLASTREADER2;

	foreach $line_first_seq(@first_sequence){
		if ( $line_first_seq=~ /^>/)    {
		$line_first_seq =~/\|(.+?)\|/;
		$gi_query = $1;
		}
	}
		
		
		
		
		my @params = ('program' => 'blastp', 'database' => "bbh_DB_$line", 'output' => "1stBBH_BLAST_$gi_query-on_DB_$line.txt", 'I' => 'T', "a" =>"30", "v"=> "1", "b"=> "0");
		my $factory = Bio::Tools::Run::StandAloneBlast->new(@params);
		my $expectvalue = 100;
		$factory->e($expectvalue);
		my $blast_report = $factory->blastall("query_seq_$i.txt");

		unless (open(BLASTREADER2, "1stBBH_BLAST_$gi_query-on_DB_$line.txt") ) {
		print "cannot open the file 1st_BLAST.txt ***\n";
		exit;
		}

		@first_blast = <BLASTREADER2>;
		close BLASTREADER2;			
	


			foreach $line1(@first_blast)  {

				if ($line1 =~ /^gi/)       {
				my $evalue= substr($line1, 76, 10);
				#print"$evalue\n\n";		
				chomp $evalue;
					if ($evalue<=1000) {
					$line1 =~/\|(.+?)\|/;
					$gi = $1;
					$line1 = $gi;                                         			
					open (SCRIVI2, ">>TABULAR_FBH_OUTPUT.xls");
					print SCRIVI2 " $line1\t";
					close SCRIVI2;

					push (@gi_first_blast, $line);
					
					#print "\n sequence $line1 passed first stage (evalue<0.005) and will be retrieved\n\n";
					#print "\n retrieving fasta...$line1\n";
					system("fastacmd -d bbh_DB_$line -s $line1 -o fasta_cmd_last_FBHit.txt");

					###########appoggio per fastacmd##########
					open(FASTAREADER2, "fasta_cmd_last_FBHit.txt");
					@fasta_cmd = <FASTAREADER2>;
					###########appoggio per fastacmd##########

					open (SCRIVI2, ">>FBH_fasta_empty_$gi_query-in-$line.txt");
					print SCRIVI2 @fasta_cmd;
					close SCRIVI2;
					}
					else{open (SCRIVI2, ">>TABULAR_BBH_OUTPUT.xls");
					print SCRIVI2 " +++ +++\t";
					close SCRIVI2;
					open (SCRIVI2, ">>TABULAR_BBH_OUTPUT.csv");
					print SCRIVI2 "Na;";
					close SCRIVI2;
					open (SCRIVI2, ">>TABULAR_FBH_OUTPUT.xls");
					print SCRIVI2 " +++ +++\t";
					close SCRIVI2;
					
					}		
				}
				else{}

		
			} 

	
		
	#print "\n>>> >>> >>> >>> >>> >>> Computing REVERSE BLAST search and SEQUENCE RETRIEVAL  *******************\n";
	my @params = ('program' => 'blastp', 'database' => "bbh_DB_query", 'output' => "Reverse_BLAST_of_query_$gi_query-in-$line.txt", 'I' => 'T', "a" => "30", "v"=> "1", "b"=> "0");
	my $factory = Bio::Tools::Run::StandAloneBlast->new(@params);
	my $expectvalue = 1000;
	$factory->e($expectvalue);
	my $blast_report = $factory->blastall("FBH_fasta_empty_$gi_query-in-$line.txt");

	unless (open(BLASTREADER2, "Reverse_BLAST_of_query_$gi_query-in-$line.txt") ) {
	print "cannot open the file 1st_BLAST.txt----\n";
	exit;
	}
	@reverse_blast = <BLASTREADER2>;
	close BLASTREADER2;
	#print @reverse_blast;

		foreach $line2(@reverse_blast){
			
			if ($line2 =~/Query= gi\|(\w+)\|ref/){
			
			$putative_orth= $1;
			chomp($putative_orth);
			#print "\n $putative_orth is a putative ortholog\n";
			}
			if ($line2 =~ /^gi/){
			$evalue2= substr($line2, 76, 10);
			#print"$evalue2\s\s";		
			chomp $evalue2;
			$line2 =~/\|(.+?)\|/;
			$gi2 = $1;
			chomp ($gi2);
			
				if (($evalue2<=0.005)&&($gi_query eq $gi2)){
				#print "\n --- --- -- -- $gi2 = $gi_query, then $putative_orth is BBH!";
				#print ".";
				push (@gi_BBHs, "$putative_orth\n");
				chomp($putative_orth);
				open (SCRIVI2, ">>TABULAR_BBH_OUTPUT.xls");
				print SCRIVI2 "$putative_orth\t";
				close SCRIVI2;
				open (SCRIVI2, ">>TABULAR_BBH_OUTPUT.csv");
				print SCRIVI2 "1;";
				close SCRIVI2;
				}
				else{
				open (SCRIVI2, ">>TABULAR_BBH_OUTPUT.xls");
				print SCRIVI2 " --- ---\t";
				close SCRIVI2;

				open (SCRIVI2, ">>TABULAR_BBH_OUTPUT.csv");
				print SCRIVI2 "Na;";
				close SCRIVI2;

				}
								

			}


		}


}
#print "\nok";



foreach $line3(@gi_BBHs){
chomp($line3);
system("fastacmd -d all_genomes_to_probe -s $line3 -o BBH_of_query.txt");
###########appoggio per fastacmd##########
open(FASTAREADER2, "BBH_of_query.txt");
@fasta_cmd2 = <FASTAREADER2>;
###########appoggio per fastacmd##########
open (SCRIVI2, ">>BBH_of_query_$gi_query.fasta");
print SCRIVI2 @fasta_cmd2;
close SCRIVI2;

}

@gi_BBHs=();	
unlink(BBH_of_query.txt);
open (SCRIVI2, ">>TABULAR_BBH_OUTPUT.xls");
print SCRIVI2 "\n";
close SCRIVI2;

open (SCRIVI2, ">>TABULAR_BBH_OUTPUT.csv");
print SCRIVI2 "\n";
close SCRIVI2;

open (SCRIVI2, ">>TABULAR_FBH_OUTPUT.xls");
print SCRIVI2 "\n";
close SCRIVI2;

@reverse_blast=();
}

print "\nCleaning up useless files (don't mind about eventual bash errors here!) \n";

system("rm *.pal");
system("rm FBH*");
system("rm *.fasta");
system("rm 1stBBH*");
system("rm Reverse*");
system("rm bbh_DB_*");
system("rm *.faa");

print "...OK\n";

print "\n\n******************          DONE           ******************\n\n\n";