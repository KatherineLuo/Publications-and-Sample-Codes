#! /usr/bin/perl -w

use strict;
use warnings;

open LIST_FILE,"find -name '*.gbk' -print|" or die ("Cannot open the find!\n");
my @file = <LIST_FILE>;
chomp(@file);
close(LIST_FILE);

my ($chromatin,$chromatin_num);
my $file_num = scalar(@file);
for(my $i=0;$i<$file_num;$i++)
{
   $chromatin = $file[$i];
   $chromatin_num = $chromatin;
   print "*****************$chromatin\n"; 
   $chromatin_num =~ s/.{4}$//i;
   open(IN,"$chromatin") or die "can't open the file $_!";
   open(GENE,">>chromatin_1.gene") or die "can't open the file!\n";
   open(CDS,">>chromatin_1.cds") or die "can't open the file!\n";
   open(EXON,">>chromatin_1.exon") or die "can't open the file!\n";
   open(GENOME,">>chromatin_1.genome") or die "can't open the file!\n";
   
   my $arr;
   my $line;
   my $line1;
   
   while(<IN>)
   {
      $line1 = $_;
      $arr .= $line1;
   }
   
   if($arr =~ /\/\/\w+/)
   {
      $arr =~ s/\/\/\w+//g;
   }
   $chromatin_num =~ s/\.\///ig;
   
   my @arr = split /\/\//,$arr;
   pop (@arr);
   #my $contig_num = scalar(@arr);
   
   ######deal with gene###############
   my $contig_mark = 0;
   foreach my $line(@arr)
   { 
      $contig_mark++;
      
      my $contig_name = "";
      my $specie = "";
      
      if($line =~ /ACCESSION/)
      {
         $line =~ /ACCESSION\s{3}(\w+)/;
	 $contig_name = $1;
      }
      #print "***$chro_name\n";
      if($line =~ /SOURCE/)
      {
         $line =~ /SOURCE\s{6,}(\w+\s+\w+)/;
         $specie = $1;
      }
      
      #//////////////deal with each gene///////////////
      while($line=~/^\s{5}gene(.*)\n(^\s{21}.*\n)*/gm)
      {
         my $ele = $&;
         my ($gene_name,$locus_tag,$stran_name,$start,$end) = gene_info($ele);
	 print GENE "$contig_name\t$gene_name\t$locus_tag\t$stran_name\t$start\t$end\n";
         
      }

      #////////////////deal with exon//////////////////////
      my $cds_num = 0;
      while($line=~/^\s{5}CDS(.*)\n(^\s{21}.*\n)*/gm)
      {
         my $ele = $&;
         my ($gene_name,$locus_tag,$stran_name,$exon_num,$note,$codon_start,$product,$pro_id,$loc) = cds($ele);
     
         my @loc = @$loc;
         print CDS "$contig_name\t$gene_name\t$locus_tag\t$stran_name\t$exon_num\t$note\t$codon_start\t$product\t$pro_id\t@loc\n";
         my $site;
         my @site;
         my $exon_mark = 0;
	 foreach my $exon(@loc)
         { 
            $exon_mark++;
	    if($exon =~ /\.\./)
	    {
	       @site = split /\.\./,$exon;
	    
	    }else{
	             $site[0] = $exon;
                     $site[1] = 0;
                  }  
            print EXON "$contig_name\t$gene_name\t$locus_tag\t$stran_name\t$pro_id\t";
	    
	    if($site[0] =~ /<|>/)
	    {
		$site[0] =~ s/<|>//g;
	    }
	    if($site[1] =~ />|</)
	    {   
		$site[1] =~ s/>|<//g;
	    }
	    print EXON "$site[0]\t$site[1]\n";
	          
         }
	 $cds_num++;
      }

      #//////////////get each contig's sequence////////////////////////////////////
      my $contig_seq = "";
      if($line =~ /ORIGIN(.*)/gs)
      {
	  $contig_seq = $1;
	  $contig_seq =~ s/\s|\d//g;
      }	 
      
      print GENOME "$specie\t$chromatin_num\t$contig_name\t$cds_num\t$contig_seq\n";
   }

   close IN;
   close GENE;
   close CDS;
   close EXON;
   close GENOME;
}

#//////////////////the main of the program has been finished/////////////////////////////

sub gene_info
{
   my ($gene) = @_;
   my ($stran_name,@gene,$gene_name,$locus_tag);
   
   $gene_name = "";
   $locus_tag = "";
   $stran_name = "";
   #$note = "";

   @gene = split /\s{21}\//,$gene;
   
   if($gene[0]=~/complement/g)
   {
      $stran_name = "crick";
   }else{
           $stran_name = "waterson";
	}
   
   my $start = "";
   my $end = "";
   
   if($gene[0] =~ /\.\./)
   {
      $gene[0] =~ /(<?\d+)\.\.(.\d+)/;
      
      $start = $1;
      $end = $2;
   }else{
           $gene[0] =~ /(\d+)/;
	   $start = $1;
	}   
   
   my (@annotation,$annotation);
   my $len = scalar(@gene);
   for (my $i=1;$i<$len;$i++)
   {
      @annotation = split /=/,$gene[$i];
      if($annotation[0] =~ /gene/)
      {
         $gene_name = $annotation[1];
         $gene_name =~ s/"|\W//g;
      }elsif($annotation[0] =~ /locus_tag/)
       {
          $locus_tag = $annotation[1];
          $locus_tag =~ s/"|\s{3,}|\n//g;
       }
        #elsif($annotation[0] =~ /note/)
        #{
	#   $note = $annotation[1];
	#   $note =~ s/"|\s{3,}|\n//gm;
	#}   
   } 
   
   return ($gene_name,$locus_tag,$stran_name,$start,$end);
}

sub cds
{
   my ($cds) = @_;
   my ($stran_name,$gene_name,$locus_tag,$note,$codon_start,$product,$pro_id);   
   $gene_name = "";
   $codon_start = "";
   $product = "";
   $pro_id = "";
   $locus_tag = "";
   $note = "";
   #$translation = "";
   
   my @cds = split /\s{21}\//,$cds;
   if($cds[0]=~/complement/g)
   {
      $stran_name = "crick";
   }else{
           $stran_name = "waterson";
	}
   $cds[0]=~s/[a-zA-Z]|\(|\)|\s{3,}//gm;

   my @loc = split /,/,$cds[0];
   my $exon_num = scalar(@loc);
   $loc[$exon_num-1] =~ s/\n//g;
   
   my $len = scalar(@cds);

   my ($annotation,@annotation);
   for (my $i=1;$i<$len;$i++)
   {
      @annotation = split /=(?!\s)/,$cds[$i];

      if($annotation[0] =~ /gene/)
      {
         $gene_name = $annotation[1];
         $gene_name =~ s/"|\W//g;
      }elsif($annotation[0] =~ /locus_tag/)
       {
          $locus_tag = $annotation[1];
          $locus_tag =~ s/"|\s{3,}|\n//gm;
       }
       elsif($annotation[0] =~ /note/)
       {
          $note = $annotation[1];
	  $note =~ s/"|\s{3,}|\n//gm;
       }elsif($annotation[0] =~ /codon_start/)
        {
           $codon_start = $annotation[1];
           $codon_start =~ s/\s|\n//g;
	}elsif($annotation[0] =~ /product/)
         {
            $product = $annotation[1];
	    $product =~ s/"|\s{3,}|\n//gm;
         }elsif($annotation[0] =~ /protein_id/)
          {
	     $pro_id = $annotation[1];
	     $pro_id =~ s/"|\s{3,}|\n//gm;
	  }#elsif($annotation[0] =~ /translation/)
	   #{
	   #   $translation = $annotation[1];
	   #   $translation =~ s/"|\s{3,}//gm;
	   #}   
   }

   return($gene_name,$locus_tag,$stran_name,$exon_num,$note,$codon_start,$product,$pro_id,\@loc);
}

#///////////////////////////////////end/////////////////////////////////////////////////
              
