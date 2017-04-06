#include <iostream.h>
#include <fstream.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include <fstream.h>
#include <iomanip.h>
#include <unistd.h>
#include <sys/types.h>

#include "defs.h"
#include "genomes.h"

/////////////////////////////////////////////////////////////
template <class T>
T max(T a, T b)
{
   return (a>b ? a : b);
}

template <class T>
T min(T a, T b)
{
   return (a>b ? b : a);
}

const char * getGeneName(const struct CDS * cds)
{
   if(!strcmp(cds->gene, ""))
      return cds->locus;
   else
      return cds->gene;
}

void freeAllOverlappingPairs(struct OverlappingPair * ptrOLPs)
{
   struct OverlappingPair * olp = ptrOLPs, * t;
   while(olp != NULL)
   {   
     //                              cout << " ********** " << endl;
      t = olp->next;	   
      if(olp->type != NULL)
	 free(olp->type);
      if(olp->codonType != NULL)
	 free(olp->codonType);

      free(olp);
      olp = t;
   }
   
   cout << "   All overlapping pairs freed...\n";
}

void walkAllOverlappingPairs(const struct OverlappingPair * ptrOLPs)
{
   const struct OverlappingPair * olp = ptrOLPs;
   int n = 0;
   while(olp != NULL)
   {  
       //cout << "   overlapping pair: " << endl;   
     printf("   %d: [ p: %s, %c, q: %s, %c, t: %s, ps: %ld, pe: %ld, qs: %ld, qe: %ld, c: %s]\n",
         ++n,
         getGeneName(olp->p),
         olp->p->strand, 
         getGeneName(olp->q),
	 olp->q->strand,
	 (olp->type == NULL ? "" : olp->type),
	 olp->ps,
	 olp->pe,
	 olp->qs,
	 olp->qe,
	 olp->codonType,
	 olp->subPCDS,
	 olp->subQCDS);
     
      olp = olp->next;
   }
}

void getPositiveSeqCodePos(char * string, int s)
{
   cout << "   ";
   for(int i=0; i<strlen(string); i++)
   {
      switch((s+i) % 3)
      {
         case 0: cout << "1"; break;
         case 1: cout << "2"; break;
         case 2: cout << "3"; break;
      }
   }   
}	

void getNegativeSeqCodePos(char * string, int s)
{
   cout << "   ";
   for(int i=0; i<strlen(string); i++)
   {
      switch((s+i) % 3)
      {
         case 0: cout << "3"; break;
         case 1: cout << "2"; break;
         case 2: cout << "1"; break;
      }
   }
}

char * revSequence(char * str)
{
   int i, l = strlen(str);
   for(i=0; i<l/2; ++i)
   {
      char c = str[l - i - 1];
      str[l - i - 1] = str[i];
      str[i] = c;
   }
   return str;
}

struct OverlappingPair * createOneOverlappingPair(struct CDS * p, struct CDS * q)
{
   struct OverlappingPair * ptrOlp = NULL;

   long maxS_S = max(p->exons->start, q->exons->start);
   long minE_E = min(p->exons->end, q->exons->end);

   if(((maxS_S < minE_E) || (maxS_S == minE_E)) && ((p->exons->start != q->exons->start) && (p->exons->end != q->exons->end)))
   {
      ptrOlp = myMalloc(struct OverlappingPair, 1);
      ptrOlp->type = NULL;
      ptrOlp->codonType = NULL;
      ptrOlp->subPCDS = NULL;
      ptrOlp->subQCDS = NULL;
      ptrOlp->next = NULL;
      
      ptrOlp->p = p;
      ptrOlp->q = q;      

      if( p->strand != q->strand)
         if(p->strand == '+')
            ptrOlp->type = saveString("-><-", 4);
         else ptrOlp->type = saveString("<-->", 4);
      else 
         ptrOlp->type = saveString("->->", 4);
      
      ptrOlp->ps = maxS_S - p->exons->start;
      ptrOlp->pe = minE_E - p->exons->start;
      ptrOlp->qs = maxS_S - q->exons->start;
      ptrOlp->qe = minE_E - q->exons->start;
      
      if(p->strand != q->strand)
      {
         int r = (ptrOlp->pe + 1) % 3;
	 if(r == 0)
         {
            int c =abs((strlen(q->CDS) - ptrOlp->qe + 1) % 3);
            switch(c)
            {
                case 0 : ptrOlp->codonType = saveString("3--3", 4); break;
                case 1 : ptrOlp->codonType = saveString("2--2", 4); break;
                case 2 : ptrOlp->codonType = saveString("1--1", 4); break;
            }
         }

	 else if(r == 1)
         {
             switch((strlen(q->CDS) - ptrOlp->qe + 1) % 3)
             {
                 case 0 : ptrOlp->codonType = saveString("2--2", 4); break;
                 case 1 : ptrOlp->codonType = saveString("1--1", 4); break;
                 case 2 : ptrOlp->codonType = saveString("3--3", 4); break;
             }
         }
   
         else if(r == 2)
         {
            switch((strlen(q->CDS) - ptrOlp->qe + 1) % 3)
            {
               case 0 : ptrOlp->codonType = saveString("1--1", 4); break;
	       case 1 : ptrOlp->codonType = saveString("3--3", 4); break;
               case 2 : ptrOlp->codonType = saveString("2--2", 4); break;
            }
         }
      }
      
      else
      {
         int c = abs((ptrOlp->pe + 1) % 3 - (ptrOlp->qe + 1) % 3); 
	 switch(c)
         {
            case 0 : ptrOlp->codonType = saveString("0--0", 4); break;
	    case 1 : ptrOlp->codonType = saveString("1--2", 4); break;
	    case 2 : ptrOlp->codonType = saveString("1--3", 4);	break;  
	 }
      }
   
      if(p->strand == '+') 
         ptrOlp->subPCDS = saveString(p->CDS + ptrOlp->ps, ptrOlp->pe - ptrOlp->ps + 1);
      else
	 ptrOlp->subPCDS = revSequence(saveString(p->CDS + strlen(p->CDS) - ptrOlp->pe - 1, ptrOlp->pe - ptrOlp->ps + 1));
      
      if(q->strand == '+') 
	 ptrOlp->subQCDS = saveString(q->CDS + ptrOlp->qs, ptrOlp->qe - ptrOlp->qs + 1);
      else
	 ptrOlp->subQCDS = revSequence(saveString(q->CDS + strlen(q->CDS) - ptrOlp->qe - 1, ptrOlp->qe - ptrOlp->qs + 1));
   }
   return ptrOlp;
}

void insertIntoOverlappingPairsTable(const struct OverlappingPair * ptrOLPs)
{
   const struct OverlappingPair * olp = ptrOLPs;
      
   while(olp != NULL)
   {  
      sprintf(sBuffer, "INSERT INTO overlappingpairs VALUES (%u, '%s', '%s', '', '%c', '%s', '', '%c', '%s', '%s', %u, %u, %u, %u, %u, '%s', '', '%s', '')",
	                   olp->p->pc->id,
                           olp->p->pc->ac, 
			   getGeneName(olp->p),
			   olp->p->strand, 
			   getGeneName(olp->q), 
			   olp->q->strand, 
                           olp->type, 
			   olp->codonType,
			   (olp->pe - olp->ps + 1),
			   olp->ps + 1, 
			   olp->pe + 1, 
			   olp->qs + 1,
			   olp->qe + 1,
                           olp->subPCDS,
			   olp->subQCDS); 

      executeSQL(sBuffer);
     
      olp = olp->next;
   }
}   

struct ContigOverlappingGenes
{
   struct Contig * c;
   struct OverlappingPair * olps;
   struct ContigOverlappingGenes * next;
};

const struct OverlappingPair * getOLP(const struct CDS * p, const struct CDS * q, const struct ContigOverlappingGenes * ptrCOLPs)
{
   const struct ContigOverlappingGenes * colp = ptrCOLPs;
   while(colp != NULL)
   {   
      const struct OverlappingPair * olp = colp->olps;
      
      while(olp != NULL)
      {
         if((olp->p == p && olp->q == q) || (olp->p == q && olp->q == p))
            return olp;

         olp = olp->next;
      }
   
      colp = colp->next; 
   }

   return NULL;
   cout << "What is the matter?\n";
}

struct OverlappingPair * getAllOverlapingPairs_modified(const struct Contig * c, FILE * outFile)
{
  struct OverlappingPair * ptrOLPs = NULL;
  int num = 0;
  if(c != NULL)
    {
      struct CDS * p = c->cds;
      struct OverlappingPair * tolp, * tolp1, * cur;
      cout << ">>> Creating all overlapping pairs on " << c->ac;

      tolp = NULL; tolp1 = NULL;
 
      while(p != NULL && p->next != NULL)
	{
	  tolp = createOneOverlappingPair(p, p->next);

	  if(tolp != NULL)
	    {
	      if(ptrOLPs == NULL)
		ptrOLPs = tolp;
	      else
		cur->next = tolp;

	      cur = tolp;
	      ++ num;
	    }
	  if(p->next->next != NULL)
	    tolp1 = createOneOverlappingPair(p, p->next->next);

	  if(tolp1 != NULL)
	    {
	      if(ptrOLPs == NULL)
		ptrOLPs = tolp1;
	      else
		cur->next = tolp1;

	      cur = tolp1;
	      ++ num;
	    }
	  p = p->next;
	}

      fprintf(outFile, "%s\t%s\t%u\t%u\t", c->ac, c->spe, c->len, c->cdsNum);
      cout << " [ " << num << " ].\n";
    }
  fprintf(outFile, "%d\t", num);

  return ptrOLPs;
}


struct OverlappingPair * getAllOverlapingPairs(const struct Contig * c, FILE * outFile)
{
   struct OverlappingPair * ptrOLPs = NULL;
   int num = 0;
   
   if(c != NULL)
   {
      struct CDS * p = c->cds;
      struct OverlappingPair * tolp, * tolp1, * cur;
      
      cout << ">>> Creating all overlapping pairs on " << c->ac;   

      while(p != NULL && p->next != NULL)
      {
	
	 tolp = createOneOverlappingPair(p, p->next);
	 
	 if(tolp != NULL)
	 {
	    if(ptrOLPs == NULL)
	       ptrOLPs = tolp;
            else
	       cur->next = tolp;

            cur = tolp;
            ++ num;
	 }
	 /*
	 insert into overlappingpairs values (581, 'NC_007714', 'SGP2_0022', '', '-', 'SGP2_0023', '', '-', '->->', '1--3', 7,  388, 395, 1, 7, '', '', '', '');
	 update overlappingpairs set pgCDS = 'ATGACTAGCACCCTGTCGCACCATGACCACGTATTCAAAAAGTTCCTCGGCGATATCGCTGTAGCTCGGGACTTTCTGGAAATCCATCTGCCGCCGCATCTGCGCAAACACTGTGATTTCAGTACTCTGGCGATGGCTTCCGGCTCGTTTATCGAAGACGATCTCAAGGGCCAGTGCTCGGGCAAGGGCGCGTGCGCTATCCATTTTTTTGCCCGTCCGTACTCCACAGCTTCTTTGAGGAGTTCAACCTCCATCGTCTTTTTCCCCAGCAACCGTTGTAGCTCTTTAATCTACTTCATGGCAGCAGTGAGTGCCGAGGCCGGAACGACCTCTTCTCCGGCGCAGACGGCGGTCAGGCTGCCATTCTGGTACTGTCGTCGCCAGGTGAACACCTGA' where contig = 'NC_007714' and pg = 'SGP2_0022' and qg = 'SGP2_0023';
	 update overlappingpairs set qgCDS = 'ATGTCGTTATTGTTACAGCGCGTCGAATGCATGAAGGAATATTCGCGATTGGCGGGTCTTGCTGAGGAACGAGAAGCGCGGGGAGAATGGAGGCAAGTTGCCGCGCTCTGGGAAAGGGCGGCGGAAGCAGGACGGCAAGTGAATCATGGCGACAAGGCCATCGCTCGGTTAGCGGCCTGCCGCCGTCGTATCGAAAATCAGGAAAACGATGACTAG' where contig = 'NC_007714' and pg = 'SGP2_0022' and qg = 'SGP2_0023';
	 
	  [ p: SGP2_0021, +, q: SGP2_0023, -, t: -><-, ps: 196, pe: 411, qs: 0, qe: 215, c: 3--3]
	  there is something wrong in genome 'NC_007714'. so when run the pro, first run this pro with mask the tolp1 part, and insert the above entry manually, then use csh to run others.
	 */
	  if(p->next->next != NULL) 
	   tolp1 = createOneOverlappingPair(p, p->next->next);

         if(tolp1 != NULL)
	   {
	     if(ptrOLPs == NULL)
               ptrOLPs = tolp1;
	     else
               cur->next = tolp1;

	     cur = tolp1;
	     ++ num;
	   }
	 
	 p = p->next;
      }
      
      fprintf(outFile, "%s\t%s\t%u\t%u\t", c->ac, c->spe, c->len, c->cdsNum); 
      
      cout << " [ " << num << " ].\n";
   }
   
   fprintf(outFile, "%d\t", num);
   
   return ptrOLPs;
}

void getOverlappingCodonType(struct OverlappingPair * olp, int& fn, int& sn, int& tn, int& an, int& nn, int& mn)
{
   if(! strcmp(olp->codonType, "1--1"))
      fn ++;
   else if(! strcmp(olp->codonType, "2--2"))
      sn ++;
   else if(! strcmp(olp->codonType, "3--3"))
      tn ++;
   else if(! strcmp(olp->codonType, "0--0"))
      an ++;
   else if(! strcmp(olp->codonType, "1--2"))
      nn ++;
   else if(! strcmp(olp->codonType, "1--3"))
      mn ++;
   
}

void getDifferlappingTypeNum(struct OverlappingPair * ptrOLPs)
{
   struct OverlappingPair * olp = ptrOLPs;

   int type[3];

   int nType [3][6];

   for(int i=0; i<3; ++i)
   {
      type[i] = 0;
      for(int j=0; j<6; ++j)
         nType[i][j] = 0;
   }  //initialize
   //cout << "Get different overlapping pairs:\n";
   char  * cTypes[] = { "-><-", "<-->", "->->"};
    int totalType[6];
    for(int i=0; i<6; ++i)
       totalType[i] = 0;
    
    for(int i=0; i<6; ++i)
    {
       for(int j=0; j<3; ++j)
          totalType[i] += nType[j][i];
	  	  
       cout << totalType[i] << "  ";
    }   
    cout << endl;
				   
}

//stat the most appearance ratio of the same amino acids
void getCodonDistribution(const char * CDS, long num[][7])
{
   for(int k=0; k<strlen(CDS)/3; ++k)
      for(int i=0; i<20; ++i)
         for(int j=0; j<aminoAcids[i].num; ++j)
            if(! strncmp(CDS + 3 * k, aminoAcids[i].codons[j], 3))
	    {
               num[i][j] ++;
	       num[i][6] ++;
            }
}

void outputSitesAndCodonPosition(const struct OverlappingPair * olp)
{
   if(olp->p->strand == '+')
      getPositiveSeqCodePos(olp->subPCDS, olp->ps); 
   else
      getNegativeSeqCodePos(olp->subPCDS, olp->ps);
   	      
   cout << olp->p->strand << endl << "   " << olp->subPCDS << endl;
   cout << "   " << olp->subQCDS << endl;
  
   if(olp->q->strand == '+')
       getPositiveSeqCodePos(olp->subQCDS, olp->qs);	    
   else
       getNegativeSeqCodePos(olp->subQCDS, olp->qs);
    
   cout << olp->q->strand << endl;
} 

char getSingleAminoAcid(const char * codon)
{
   for(int i=0; i<AMINO_ACID_NUM; ++i)
   {
      for(int j=0; j<aminoAcids[i].num; ++j)
         if(!strncmp(codon, aminoAcids[i].codons[j], 3))
            return aminoAcids[i].aa;
   }

  return '-';
}

char * getOlpProteinSequence(char * CDS)
{
   int i, l = strlen(CDS);
   char * AAs = myMalloc(char, l / 3 + 1); // included '*'
   
   for(i=0; i<l/3; ++i)
      AAs[i] = getSingleAminoAcid(CDS + 3 * i);

   AAs[i] = '\0';
   return AAs;
}
//supply lacking codons 
void getOverlapCodonDistribution(const struct OverlappingPair * olp, long m[][7], long n[][7])
{
   char * subPCDS;
   char * subQCDS;
   
   char * olpPAAs;
   char * olpQAAs;
  
   char * olpPAAsSelf;
   char * olpQAAsSelf;
   
   //cout << "   " << olp->type << endl; 
   if(olp->p->strand == '+')
   {
      int delta = olp->ps % 3;
      subPCDS = saveString(olp->p->CDS + (olp->ps - delta), olp->pe - olp->ps + 1 + delta);
      olpPAAs = saveString(olp->p->AAs + ((olp->ps - delta) / 3), (olp->pe + delta)/ 3 - olp->ps / 3);
   }
   else // olp->p->strand = '-' 
   {
      int delta = olp->ps % 3;
      subPCDS = saveString(olp->p->CDS + (strlen(olp->p->CDS) - olp->pe - 1 ), olp->pe - olp->ps + 1 + delta);
      olpPAAs = saveString(olp->p->AAs + strlen(olp->p->AAs) - (olp-> pe / 3), (olp->pe + delta)/ 3 - (olp->ps / 3));
   }
   
   olpPAAsSelf = getOlpProteinSequence(subPCDS);
   /*cout << " p:" << getGeneName(olp->p) 
        << ": " << olp->ps << ", " << olp->pe << ", p CDS length: " << strlen(olp->p->CDS) 
	<< ", strand: " << olp->p->strand << endl
        << "   " << olp->subPCDS << "*****" << endl 
	<< "   " << (olp->p->strand == '+' ? subPCDS : revSequence(subPCDS)) << endl; 
   
   cout << "   " << olpPAAs << endl << "   " << olpPAAsSelf<< endl; 
   */
   if(olp->q->strand == '+')
   {
      int delta = 3 - (olp->qe + 1) % 3;
      subQCDS = saveString(olp->q->CDS + (olp->qs), olp->qe - olp->qs + 1 + delta);
      olpQAAs = saveString(olp->q->AAs + (olp->qs / 3), (olp->qe + delta) / 3 - (olp->qs / 3) + 1);
      
   } 
   else // olp->q->starnd = '-'
   {
      int delta = 3 - (olp->qe + 1) % 3;
      subQCDS = saveString(olp->q->CDS + (strlen(olp->q->CDS) - olp->qe - 1 - delta), olp->qe + delta - olp->qs + 1 );
      olpQAAs = saveString(olp->q->AAs + strlen(olp->q->AAs) - (olp->qe / 3) - 1, (olp->qe + delta) / 3 - (olp->qs / 3) + 1);
      
   } 
   
   olpQAAsSelf = getOlpProteinSequence(subQCDS);
   /*cout << " q:" << getGeneName(olp->q) 
        << ": " << olp->qs << ", " << olp->qe << ", q CDS length: " << strlen(olp->q->CDS) 
	<< ", strand: " << olp->q->strand << endl 
        << "   " << subQCDS << "*****" <<endl 
	<< "   " << (olp->q->strand == '+' ? subQCDS : revSequence(subQCDS)) << endl; 

   cout << "   " << olpQAAs << endl << "   " << olpQAAsSelf << endl << endl;
   */
   sprintf(sBuffer, "UPDATE overlappingpairs SET polpAAs='%s', qolpAAs='%s'  WHERE pg='%s' AND qg='%s'", olpPAAs, olpQAAs, getGeneName(olp->p), getGeneName(olp->q));
   executeSQL(sBuffer);
    
   getCodonDistribution(subPCDS, m);
   getCodonDistribution(subQCDS, m);
   
   getCodonDistribution(olp->p->CDS, n);
   getCodonDistribution(olp->q->CDS, n);
   
   free(subPCDS);
   free(subQCDS); 
   free(olpPAAs);
   free(olpQAAs);
   free(olpPAAsSelf);
   free(olpQAAsSelf);
   
}	

void outputAnAminoAcidCodonDistribution(int i, long * num, FILE * outFile)
{
   fprintf(outFile, "%c [%6u] : ", aminoAcids[i].aa, num[6]);
   for(int j=0; j<aminoAcids[i].num; j++)
   {
      if(j)
         fprintf(outFile, ", ");   	 
     
      if(num[6] == 0)     
         fprintf(outFile, "%3s - %6u - %6.1f%%", aminoAcids[i].codons[j], num[j], 0.0);	 
      else 
         fprintf(outFile, "%3s - %6u - %6.1f%%", aminoAcids[i].codons[j], num[j], float(num[j]) / num[6] * 100.0 );	 
   }
   fprintf(outFile, "\n");
}

float chiSquareTestValue(int i, long * n, long * m)
{
   float statistic1 = 0.0;
   float statistic2 = 0.0;
   float frequency1 = 0.0;
   float frequency2 = 0.0;
   
   if(aminoAcids[i].num > 1)
      for(int j=0; j<aminoAcids[i].num; j++)
      {
	 frequency1 = float(n[6] * (n[j] + m[j])) / (n[6] + m[6]);
	 frequency2 = float(m[6] * (n[j] + m[j])) / (n[6] + m[6]);
         
	 statistic1 += (n[j] - frequency1) * (n[j] - frequency1) / frequency1;
         statistic2 += (m[j] - frequency2) * (m[j] - frequency2) / frequency2; 
      }
   
   return statistic1 + statistic2;
   
}

void getThresholdValues(int i, float * tv01, float * tv05)
{
   int   df = aminoAcids[i].num - 2; //df>=1
   float zeroFive[5] = {3.84, 5.99, 7.81,  9.49,  11.07};
   float zeroOne[5]  = {6.63, 9.21, 11.34, 13.28, 15.09};
   
   * tv01 = zeroOne[df];
   * tv05 = zeroFive[df];
}
	
////////////////////////////////////////////////////////////////////
// Wilcoxon signed-ranks test for the median difference
float getWilcoxonTestValue(int i, long * m, long * n, float * tv01, float * tv05)
{
   float x[6];
   int k = 0, signs[6];
   for(int j=0; j<aminoAcids[i].num; ++j)
   {
      int t = m[j] - n[j];
      if(t>0)
      {
         x[k] = t;
	 signs[k++] = 1;
      }
      else if(t<0)
      {
         x[k] = -t;
         signs[k++] = -1;
      }	 
   }

   // simple sorting method due to too small numbers
   for(int u=0; u<k; ++u)
   {
      for(int v=u+1; v<k; ++v)
      {
         if(x[v] < x[u])
	 {
            int tmp = int(x[u]);
            x[u] = x[v];
            x[v] = tmp;
	 }
      }
   }
   
   // get positive sum...
   for(int u=0; u<k; ++u)
   {
      float sum = u + 1;
      int v;      
      for(v=u+1; v<k; ++v)
      {
         if(x[v] > x[u])
            break;
	 else
	    sum += v + 1;	 
      }

      for(int w=u; w<v-1; ++v)
         x[w] = sum / (v - u);
   }

   float sum = 0;
   for(int w=0; w<k; ++w)
      if(signs[w] == 1) // positive
         sum += x[w];

   float upperTailValues01[] = {};

   return sum;
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////

void outputCodonDistribution(const char * ac, const char * type, long m[][7], long n[][7], FILE * outFile)
{
   float tv01, tv05;	
   for(int i=0; i<20; i++)
   { 
      if(aminoAcids[i].num <= 1)
         continue;	      
      // outputAnAminoAcidCodonDistribution(i, m[i], outFile);
      // outputAnAminoAcidCodonDistribution(i, n[i], outFile);
      getThresholdValues(i, &tv01, &tv05);
      float statistics = chiSquareTestValue(i, m[i], n[i]); 
      fprintf(outFile, "\"%s\"\t\"%s\"\t\"%c\"\t\"%c\"\t%f\t%u\t%f\t%f\n",
		       ac, type,
                       aminoAcids[i].aa,
                       (statistics > tv05 ? 'Y' : 'N'),
		       statistics,
		       aminoAcids[i].num - 1,
		       tv05, tv01);
   }

   //fprintf(outFile, "\n");
}
            
void setInitialValues(long num[][7])
{
   for(int i=0; i<20; ++i)
      for(int j=0; j<7; ++j)
         num[i][j] = 0;
}

long getAllSitesAndCodonPosition(const struct Contig * c, struct OverlappingPair * ptrOLPs, FILE * outFile, FILE * outFile2)
{
   long hhNum[20][7], hhNumTotal[20][7];
   long ttNum[20][7], ttNumTotal[20][7];
   long coNum[20][7], coNumTotal[20][7];
 
   setInitialValues(hhNum);
   setInitialValues(ttNum);
   setInitialValues(coNum);
   setInitialValues(hhNumTotal);
   setInitialValues(ttNumTotal);
   setInitialValues(coNumTotal);
    
   struct OverlappingPair * olp = ptrOLPs;
   long codingSize = 0;
   long overlapSize = 0;
   
   while(olp != NULL)
   {
      overlapSize += strlen(olp->subPCDS) + strlen(olp->subQCDS); // just each overlapping pair subPCDS and subQCDS
      
      codingSize += strlen(olp->p->CDS) + strlen(olp->q->CDS);
      
      if(!strcmp(olp->type, "-><-"))
	 getOverlapCodonDistribution(olp, hhNum, hhNumTotal); // after suppplying the lacked amino acids and translation
      else if(!strcmp(olp->type, "<-->"))
	 getOverlapCodonDistribution(olp, ttNum, ttNumTotal);
      else if(!strcmp(olp->type, "->->"))
         getOverlapCodonDistribution(olp, coNum, coNumTotal);     
      
      olp = olp->next;
   } 
   
   outputCodonDistribution(c->ac, "-><-", hhNum, hhNumTotal, outFile); //after chiSquare examing, each amino acid distribution
   outputCodonDistribution(c->ac, "<-->", ttNum, ttNumTotal, outFile);
   outputCodonDistribution(c->ac, "->->", coNum, coNumTotal, outFile);   
   
   fprintf(outFile2, "%u\t%u", codingSize, overlapSize);
   return overlapSize; // total overlapping pair size
}   

void outputInFasta(const struct CDS * cds, FILE * outFile, int len = 60)
{
      //cout <<  cds->pc->ac << ", " << getGeneName(cds) << ", " <<  cds->start << ", " << strlen(cds->AAs) - 1 << endl;

      fprintf(outFile, ">%s|%s|%lu|%lu\n", 
                     cds->pc->ac,
	             getGeneName(cds),
		     cds->start,
		     strlen(cds->AAs) - 1); 
   
      int l = strlen(cds->AAs);   
      int n = 0;
   
      while((n<l-1) && (l > 0)) 
      {
         fputc(cds->AAs[n++], outFile);
      
         if(n % len == 0 || n == l-1)
            fprintf(outFile, "\n");
      } //each line just contains 60 words
}   
     

struct QSBlastRelation
{
   char * c;
   char * g;
   int len;
   int start;
};

struct QSBlastRelation * analyzeStr(char * str)
{
   struct QSBlastRelation * qs = myMalloc(struct QSBlastRelation, 1);
   char * p = str;
   char * q = strstr(str, "|"); //the first position that "|" emerges in string p.
   qs->c = saveString(p, q - p);

   p = q + 1;
   q = strstr(p, "|");
   
   qs->g = saveString(p, q - p);
   
   p = q + 1;
   sscanf(p, "%u|%u", &(qs->start), &(qs->len));
  
   return qs;
}

void freeQSBlastRelation(struct QSBlastRelation * qs)
{
   if(qs->c != NULL)
      free(qs->c);
   
   if(qs->g != NULL)
      free(qs->g);

   free(qs);
}

//void insertIntoBlastTable(char * str, FILE * outFile)
void outputIntoFile(char * str, FILE * outFile)
{
   char qstr[80];
   char sstr[80];
   
   float identity;
   int al;
   int mis, gap, qs, qe, ss, se;
   char evalue[40];
   int score;

   sscanf(str, "%s %s %f %u %u %u %u %u %u %u %s %u", 
	        qstr, sstr, &identity, &al, &mis, &gap, &qs, &qe, &ss, &se, evalue, &score);
   
   struct QSBlastRelation * qQS = analyzeStr(qstr);
   struct QSBlastRelation * sQS = analyzeStr(sstr);
   
   if(strcmp(qQS->c, sQS->c) || qQS->start != sQS->start)  //isn't itself
   {
       //sprintf(sBuffer, "INSERT INTO blast VALUES ('%s', '%s', %u, %u, '%s', '%s', %u, %u, %f, %u, %u, %u, %u, %u, %u, %u, %s, %u)",
       //	    qQS->c, qQS->g, qQS->len, qQS->start, sQS->c, sQS->g, sQS->len, sQS->start,
       //    identity, al, mis, gap, qs, qe, ss, se, evalue, score);
       //executeSQL(sBuffer);

       fprintf(outFile, "%s\t%s\t%u\t%u\t%s\t%s\t%u\t%u\t%f\t%u\t%u\t%u\t%u\t%u\t%u\t%u\t%s\t%u", 
	       qQS->c, qQS->g, qQS->len, qQS->start, sQS->c, sQS->g, sQS->len, sQS->start, identity, al, mis, gap, qs, qe, ss, se, evalue, score);
       fprintf(outFile, "\n");
   }
   
   freeQSBlastRelation(qQS);
   freeQSBlastRelation(sQS);
}

//The whole file:Escherichia_coli.fasta(pid, pl) blast against the database:fasta.out
void insertIntoBlastTwoDBTable(char * str, FILE * outFile)
{
   char sstr[80];
   
   float identity;
   int qpid, qpl, al;
   int mis, gap, qs, qe, ss, se;
   char evalue[40];
   int score;

   sscanf(str, "%u|%u %s %f %u %u %u %u %u %u %u %s %u", 
	        &qpid, &qpl, sstr, &identity, &al, &mis, &gap, &qs, &qe, &ss, &se, evalue, &score);
   
   struct QSBlastRelation * sQS = analyzeStr(sstr);
   
   sprintf(sBuffer, "INSERT INTO blastTwoDB VALUES (%u, %u, '%s', '%s', %u, %u, %f, %u, %u, %u, %u, %u, %u, %u, %s, %u)",
   	                         qpid, qpl, sQS->c, sQS->g, sQS->len, sQS->start, identity, al, mis, gap, qs, qe, ss, se, evalue, score);
   executeSQL(sBuffer);

   //fprintf(outFile, "%u\t%u\t'%s'\t'%s'\t%u\t%u\t%f\t%u\t%u\t%u\t%u\t%u\t%u\t%u\t%s\t%u\n", 
   //	                         qpid, qpl, sQS->c, sQS->g, sQS->len, sQS->start, identity, al, mis, gap, qs, qe, ss, se, evalue, score);
   freeQSBlastRelation(sQS);

}   

void deleteFromBlastTwoDB()
{
   sprintf(sBuffer, "SELECT A.sc, A.sg FROM blastTwoDB AS A, overlappingpairs AS B WHERE A.sc=B.contig AND (A.sg=B.pg OR A.sg=B.qg)");
   executeSQL(sBuffer);
   mysqlResult = mysql_store_result(mysqlConnection);
   while((mysqlRow = mysql_fetch_row(mysqlResult)) != NULL)
   {
      sprintf(sBuffer, "DELETE FROM blastTwoDB WHERE sc='%s' AND sg='%s'", mysqlRow[0], mysqlRow[1]);
      executeSQL(sBuffer);
   }
   mysql_free_result(mysqlResult);
}   

void blastSearchDB(const char * prog, const char * db, const struct CDS * cds, float exp = 0.00001)
{
   pid_t procID = getpid();
   sprintf(sBuffer, "/tmp/blast%u.seq", procID);
   
   FILE * fp = fopen(sBuffer, "wt");
   outputInFasta(cds, fp); //format sequencre
   fclose(fp); 
		   
   sprintf(sBuffer, "blastall -p %s -i /tmp/blast%u.seq -d %s -e %f -m 8 > /tmp/blast%u.out",
		    prog, procID, db, exp, procID);
   system(sBuffer); // do blast search
   
   sprintf(sBuffer, "/tmp/blast%u.out", procID);
   fp = fopen(sBuffer, "rt");
   
   //cout << ">>> Query seq: " << cds->pc->ac << "|" << getGeneName(cds) << " [" << strlen(cds->AAs) << "]\n";

   while(!feof(fp))
   {
       //if(fgets(sBuffer, 255, fp) != NULL)
       //insertIntoBlastTable(sBuffer);
      //cout << sBuffer;
   }
   fclose(fp); 
}

void doBlastSearch(const struct Contig * ptrContigs)
{
   //executeSQL("DELETE FROM blast");   

   const struct Contig * c = ptrContigs;
  
   int num = 0;
   
   while(c != NULL)
   {
      const struct CDS * cds = c->cds;
      
      while(cds != NULL)
      {  	      
	 blastSearchDB("blastp", "db/B.ppb", cds);
         
            if(++num > 10) // for test
	      break;
	 
         cds = cds->next;   
      }	 
      
      c = c->next;
   }
}

void blastSearchInTwoDB(const char * ac, const char * prog, const char * db, float exp = 0.0001)
{
   pid_t procID = getpid();
   sprintf(sBuffer, "/dahome/olp/luoyq/bin/temp_data/blast-2.2.14/bin/blastall -p %s -a 3 -i fastaDB/%s.fasta -d %s -e %f -m 8 > /tmp/blast%u.out", prog, ac, db, exp, procID);
   system(sBuffer); // do blast search
   cout << sBuffer << " ---- " << endl;   
   sprintf(sBuffer, "/tmp/blast%u.out", procID);
   
   cout << "the proID is: " << procID << endl;

   //sprintf(sBuffer, "/tmp/blast17455.out");
   
   char * name = myMalloc(char, 100);
   FILE * fp = fopen(sBuffer, "rt");
   sprintf(name, "/dahome/olp/luoyq/bin/blast/%s.blast", ac);
   cout << "name:  ---- " << name << " ---- \n";
   cout << "ac:  ---- " << ac << " ---- \n";
   FILE * outFile = fopen(name, "wt");

   int i = 0;
   while(!feof(fp))
   {
      if(fgets(sBuffer, 255, fp) != NULL)	
	//insertIntoBlastTable(sBuffer);
	outputIntoFile(sBuffer, outFile);
   }
   free(name);
   fclose(fp);   
   fclose(outFile);
}

///////////////////////////////////////////////////////////////////

struct Homology
{
   struct CDS * q;	
   float identity, e, s, ap;
   int al, pl, ql; //pl:alignment length of query gene; ql: alignment length of subject gene
   struct Homology * next;
};
	
const struct CDS * getCDS(const Contig * pcs, char * ac, char * gene, int start)
{
   const struct Contig * c = pcs;
   while(c != NULL)
   {
      //cout << "*****************, " << c->ac << ", " << ac << endl;
      if(!strcmp(c->ac, ac))
      {	      
         const struct CDS * cds = c->cds;
         while(cds != NULL)
         {
	     if(cds->start == start)
		 return cds;
	   
	     cds = cds->next;
	 }	     
      } 
      
      c = c->next;
   }
   //cout << "################, " << ac << ", " << c->ac << endl;
   return NULL;
}


char * acSet()
{
  MYSQL_RES * ret;
  MYSQL_ROW row;
  
  sprintf(sBuffer, "SELECT contig FROM olp_species where sp_mark = 1");
  //sprintf(sBuffer, "SELECT contig FROM olp_species where taxo_gbk LIKE \"%%Proteobacteria%%\" AND sp_mark = 1");
  //sprintf(sBuffer, "SELECT contig FROM olp_species where taxonomy not like \"%%Proteobacteria%%\" and taxonomy not like \"%%Firmicutes%%\" and taxonomy != ''");
  //sprintf(sBuffer, "SELECT contig FROM olp_species where taxo_gbk like \"%%Gamma%%\" AND sp_mark = 1");

  executeSQL(sBuffer);
  ret = mysql_store_result(mysqlConnection);

  int num = mysql_num_rows(ret);
  int len = num * 20;
  char * acset = myMalloc(char, len);
  strcpy(acset, "");

  for(int i=0; i < num && (row = mysql_fetch_row(ret)) != NULL; i++)
    {
      if(i == 0)
        sprintf(acset, "(\'%s\', ", row[0]);
      else if(i == num - 1)
        sprintf(acset, "%s\'%s\')", acset, row[0]);
      else
        sprintf(acset, "%s\'%s\', ", acset, row[0]);
    }

  mysql_free_result(ret);

  return acset;
}

// using bidirectional comparison to find homolgies
struct Homology * getHomology(const struct CDS * p, const struct Contig * pcs, float iden, float evl)
{
  const char * gs = acSet();

    sprintf(sBuffer, "SELECT identity, al, qe-qs+1, se-ss+1, evalue, sc, sg, sgs, sgl, score FROM blast WHERE qc='%s' AND qgs=%u AND identity>%f AND evalue<%f AND sc IN %s GROUP BY sc ORDER BY evalue ASC", p->pc->ac, p->start, iden, evl, gs);
   executeSQL(sBuffer);
   mysqlResult = mysql_store_result(mysqlConnection);
   
   struct Homology * ptrHomology = NULL, * pho, * ho;
   while((mysqlRow = mysql_fetch_row(mysqlResult)) != NULL)
   {
      MYSQL_RES * ret;
      MYSQL_ROW row;      
      sprintf(sBuffer, "SELECT identity, al, qe-qs+1, se-ss+1, evalue, sc, sg, sgs, score FROM blast WHERE qc='%s' AND qgs=%u AND identity>%f AND evalue<%f AND sc IN %s GROUP BY sc ORDER BY evalue ASC", mysqlRow[5], atoi(mysqlRow[7]), iden, evl, gs);
      executeSQL(sBuffer);
      ret = mysql_store_result(mysqlConnection);
      while((row = mysql_fetch_row(ret)) != NULL)
      {
	 if(!strcmp(row[5], p->pc->ac) && atoi(row[7]) == p->start)
	 {
            ho = myMalloc(struct Homology, 1);
            ho->q = (struct CDS *) getCDS(pcs, mysqlRow[5], mysqlRow[6], atoi(mysqlRow[7]));
	    
	    if(ho->q == NULL)
	    {
	      //cout <<  "Ho.q = NULL\n";
               break;
	    }

            ho->identity = (atoi(mysqlRow[0]) + atoi(row[0])) / 2;
            ho->al = (atoi(mysqlRow[1]) + atoi(row[1])) / 2;
            ho->pl = (atoi(mysqlRow[2]) + atoi(row[2])) / 2;
            ho->ql = (atoi(mysqlRow[3]) + atoi(row[3])) / 2;
            ho->e  = (atof(mysqlRow[4]) + atof(row[4])) / 2;
	    ho->s = (atof(mysqlRow[9]) + atof(row[8])) / 2;
	    ho->ap = (float (ho->al + 1) / strlen(p->AAs));
            ho->next = NULL;
 
	    if(ptrHomology == NULL)
               ptrHomology = ho;
            else
	       pho->next = ho;
            pho = ho;

	 }
      }
      mysql_free_result(ret);
   }     
   
   mysql_free_result(mysqlResult);
   
   return ptrHomology;
}

int analyzeHomologousPair(const struct OverlappingPair * olp, const struct Homology * hps, const struct Homology * hqs, const struct ContigOverlappingGenes * ptrCOLPs)
{
   const struct Homology * hp = hps;
   int num = 0;
   
   while(hp != NULL)
   {
      const struct Homology * hq = hqs;
      while(hq != NULL)
      {
         if(!strcmp(hp->q->pc->ac, hq->q->pc->ac))
	 {
	     const struct OverlappingPair * holp = getOLP(hp->q, hq->q, ptrCOLPs);
	     sprintf(sBuffer, "INSERT INTO olp_comparison_all VALUES (0, '%s', '%s', '%c', '%s', '%c', %lu, %lu, %u, %u, '%s', '%s', %u",
			   olp->p->pc->ac, getGeneName(olp->p), olp->p->strand, getGeneName(olp->q), olp->q->strand, 
			   olp->p->start, olp->q->start, strlen(olp->p->CDS), strlen(olp->q->CDS),
			   olp->type, olp->codonType, olp->pe - olp->ps + 1);
	     //cout << sBuffer << endl;
             
	     if(holp != NULL)	     
	     {
	      	 ++ num; 
		 
		 sprintf(sBuffer, "%s, 0, '%s', '%s', '%c', '%s', '%c', %lu, %lu, %u, %u, '%s', '%s', %u, %g, %f, %f, %f)", 
				  sBuffer, hp->q->pc->ac, getGeneName(hp->q), hp->q->strand, 
				  getGeneName(hq->q), hq->q->strand,
				  hp->q->start, hq->q->start, strlen(hp->q->CDS), strlen(hq->q->CDS),
			          holp->type, holp->codonType, holp->pe - holp->ps + 1, (hp->e + hq->e) / 2,
			          (hp->s + hq->s) / 2, min(hp->ap, hq->ap), (hp->identity + hq->identity) / 2);
		 executeSQL(sBuffer);

	     }
             else
	     {
		 sprintf(sBuffer, "%s, 0, '%s', '%s', '%c', '%s', '%c', %lu, %lu, %u, %u, '', '', %u, %g, %f, %f, %f)", 
				  sBuffer, hp->q->pc->ac, getGeneName(hp->q), hp->q->strand, 
				  getGeneName(hq->q), hq->q->strand,
				  hp->q->start, hq->q->start, strlen(hp->q->CDS), strlen(hq->q->CDS), 
          			  ((hq->q->start > hp->q->start) ? abs(hp->q->start - hq->q->start - strlen(hp->q->CDS) + 1) : abs(hp->q->start - hq->q->start - strlen(hq->q->CDS) + 1)), (hp->e + hq->e) / 2,
		     	          (hp->s + hq->s) / 2, min(hp->ap, hq->ap), (hp->identity + hq->identity ) / 2);
		 executeSQL(sBuffer);

	     }
	 }
         hq = hq->next;	 
      }

      hp = hp->next;
   }

   return num;
}

void freeAllHomologies(struct Homology * ptrHs)
{
   struct Homology * h = ptrHs;
   while(h != NULL)
   {
      struct Homology * t = h->next;
      free(h);
      h = t;      
   }
}

void getOverlappingHomology(const struct Contig * pcs, const struct ContigOverlappingGenes * ptrCOLPs, float iden=0.40, float evl=0.0001)
{
   int num = 0;
   int nZero = 0;
   int nOne = 0;
   int nTwo = 0;
   int nOverlapped = 0;
   
   //executeSQL("DELETE FROM olp_comparison_all");
   cout << "\n>>> Comparing olps...\n"; 
   const struct ContigOverlappingGenes * colp = ptrCOLPs;
   while(colp != NULL)
   {   
      cout << "    Contig: " << colp->c->ac << "...\n";
      
      const struct OverlappingPair * olp = colp->olps;
      int i = 0;
      while(olp != NULL)
      {
	 struct Homology * hps = getHomology(olp->p, pcs, iden, evl);
         struct Homology * hqs = getHomology(olp->q, pcs, iden, evl);
         
	 if(hps == NULL && hqs == NULL)
	 {
	     sprintf(sBuffer, "INSERT INTO olp_comparison_all  VALUES (0, '%s', '%s', '%c', '%s', '%c', %lu, %lu, %u, %u, '%s', '%s', %u",
			   olp->p->pc->ac, getGeneName(olp->p), olp->p->strand,
			   getGeneName(olp->q), olp->q->strand,
			   olp->p->start, olp->q->start, strlen(olp->p->CDS), strlen(olp->q->CDS),
			   olp->type, olp->codonType, olp->pe - olp->ps + 1);
	     sprintf(sBuffer, "%s, 0, '', '', '', '', '', 0, 0, 0, 0, '', '', 0, 0, 0, 0, 0)", sBuffer);

	     executeSQL(sBuffer);
	     	     
	     ++ nZero;
	 }
         else if(hps == NULL || hqs == NULL)
	 {
	    ++ nOne;
	    
	    const struct Homology * h;
	    int isP = 1;
	    if(hps == NULL)
	    {
	       isP = 0;
   	       h = hqs;
	    }
	    else
	       h = hps;

            while(h != NULL)
	    {
                sprintf(sBuffer, "INSERT INTO olp_comparison_all  VALUES (0, '%s', '%s', '%c', '%s', '%c', %lu, %lu, %u, %u, '%s', '%s', %u",
		           olp->p->pc->ac, getGeneName(olp->p), olp->p->strand,
			   getGeneName(olp->q), olp->q->strand,
			   olp->p->start, olp->q->start,  strlen(olp->p->CDS), strlen(olp->q->CDS),
			   olp->type, olp->codonType, olp->pe - olp->ps + 1);
	       if(isP)
		 sprintf(sBuffer, "%s, 0, '%s', '%s', '%c', '', '', %lu, 0, %u, 0, '', '', 0, %g, %f, %f, %f)",
			  sBuffer, h->q->pc->ac, getGeneName(h->q), h->q->strand, h->q->start, strlen(h->q->CDS), h->e, h->s,  h->ap, h->identity); //
	       else
		 sprintf(sBuffer, "%s, 0, '%s', '', '', '%s', '%c', 0, %lu, 0, %u, '', '',  0, %g, %f, %f, %f)",
			 sBuffer, h->q->pc->ac, getGeneName(h->q), h->q->strand, h->q->start, strlen(h->q->CDS), h->e, h->s, h->ap, h->identity);
	       executeSQL(sBuffer);
   	       
	       h = h->next;
	    }
	    	    
	 }
         else
	 {
            ++ nTwo;
	    nOverlapped += analyzeHomologousPair(olp, hps, hqs, ptrCOLPs);
	 }
      
         freeAllHomologies(hps);
         freeAllHomologies(hqs);
      
         olp = olp->next;
      }

      colp = colp->next;
   }
   
   printf(">>> nZero = %u, nOne = %u, nTwo = %u, nOverlapped = %u\n", nZero, nOne, nTwo, nOverlapped);
}

void createAllOverlappingPairs(const struct Contig * ptrContigs, const char * ac)
{
  time_t start, end;
  time(&start);

    //doBlastSearch(ptrContigs);
    //return;
    //blastSearchInTwoDB(ac, "blastp", "temp_data/bacteria");
    //deleteFromBlastTwoDB(); //delete genes that are not overlapped.
    //return;
   
   FILE * outFile = fopen("bacteria_codonDistribution_3.out", "wt");
   fprintf(outFile, "Contig\tType\tAA\tSign\tStatistics\tdf\tCV05\tCV01\n");
   
   FILE * outFile2 = fopen("bacteria_statistics_3.out", "wt");
   fprintf(outFile2, "Contig\tGenome\tSize\tcdsNum\toverlapNum\tcodingSize\toverlapSize\tolpPercent\ttotalcodingSize\ttotalPercent\n");
  
   //FILE * outFile3 = fopen("temp_data/13gamma.fasta", "wt");
           
   struct ContigOverlappingGenes * ptrContigCOLPs = NULL;
   struct ContigOverlappingGenes * pcolp, * colp;
   
   const struct Contig * c = ptrContigs;
   
   while(c != NULL)
   {
      long cdsLength = 0;
      const struct CDS * cds = c->cds;
      
      while(cds != NULL)
      {  
	//outputInFasta(cds, outFile3); // ouput AAS in FASTA format
	 
	 cdsLength += strlen(cds->CDS);
	 
         cds = cds->next;   
      }	 

      colp = myMalloc(struct ContigOverlappingGenes, 1);
      colp->c = (struct Contig *) c;
      
      cout << "------- " << colp->c->ac << endl;
      colp->olps = getAllOverlapingPairs_modified(c, outFile2);

      colp->next = NULL;  
      
      //walkAllOverlappingPairs(colp->olps);
      //insertIntoOverlappingPairsTable(colp->olps);
      
      //long overlapSize = getAllSitesAndCodonPosition(c, colp->olps, outFile, outFile2);
      //double overlapPercent = double (overlapSize) / c->len * 100;
      //double codingPercent = double (cdsLength) / c->len * 100;
      //fprintf(outFile2, "(%f%%)\t%u(%f%%)\n", overlapPercent, cdsLength, codingPercent);
      //fprintf(outFile2, "\t%f%\t%u\t%f%\n", overlapPercent, cdsLength, codingPercent);

      //getDifferlappingTypeNum(colp->olps);
      //outputSitesAndCodonPositionsition(colp->olps);

      if(ptrContigCOLPs == NULL)
	 ptrContigCOLPs = colp;
      else
	 pcolp->next = colp;

      pcolp = colp;      

      c = c->next;
   }

   getOverlappingHomology(ptrContigs, ptrContigCOLPs);

   // free objects
   pcolp = ptrContigCOLPs;
   while(pcolp != NULL)
   {
      colp = pcolp->next;	   
      freeAllOverlappingPairs(pcolp->olps);
      free(pcolp);
      pcolp = colp;
   }

   fclose(outFile);
   fclose(outFile2);
   //fclose(outFile3);

   time(&end);
   printf("The start time:\t\t%s", ctime(&start));
   printf("The end time:\t\t%s", ctime(&end));
   printf("The program lasted:\t%ld\nc", end - start);

}

//////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////
/******************************************************************************************************
CREATE TABLE blast (
   qc varchar(30) not null, //quired contig
   qg varchar(50), //quired gene
   qgl int(6), //quired gene length
   qgs int(10), //start
   sc varchar(30) not null, //subject contig
   sg varchar(50), //subject gene
   sgl int(6), //length
   sgs int(10), //start
   identity float (5), //identity between qg and sg
   al int (6), //alignment length
   mis int (5) , //missing length
   gap int (5), 
   qs int (10) , //quired gene alignment start 
   qe int (10) , //quired gene alignment end
   ss int (10),
   se int (10),
   evalue float(10, 8), //blast value
   score int(10), //blast score
   index qc(qc),
   index qcg(qc, qg),
   index qcs(qc, qg, qgs)) TYPE=MyISAM;

CREATE TABLE blastTwoDB (
   qpid int(6) not null,
   qpl int(10),
   sc varchar(30) not null, 
   sg varchar(50), 
   sgl int(6),
   sgs int(10),
   iden double, 
   al int (6), 
   mis int (5),
   gap int (5), 
   qs int (10),  
   qe int (10),
   ss int (10),
   se int (10),
   evalue double, 
   score int(10),
   index qpid(qpid),
   index qcs(sc, sg)) TYPE=MyISAM; 

CREATE TABLE searchTwoDB (
   contig varchar(10) not null,
   gene varchar(10) not null, 
   auto_pfamseq varchar(50) not null,
   index cg(contig, gene),
   index au(auto_pfamseq)) TYPE=MyISAM;

CREATE TABLE olp_comparison_all (
   qc  varchar(30) not null, //quired contig
   qpg varchar(20), //quired contig p olp gene 
   qpgs varchar(5), //strand
   qqg varchar(20), //quired contig q olp gene
   qqgs varchar(5), 
   qgs int(10), //olp p start
   qqs int(10), //olp q start
   qpl int(6), //olp p CDS length
   qql int(6), //olp q CDS length
   qt  char(4), //olp type
   qct char(4), //olp codon type
   ql  int(10), //length of olp
   sc  varchar(30) not null, //subject contig
   spg varchar(20), //subject contig p gene
   spgs varchar(5), //strand
   sqg varchar(20), //subject contig q gene
   sqgs varchar(5), 
   sps int(10), //homology olp p start
   sgs int(10), //homology olp q start
   spl int(6), //subject olp p length 
   sql int(6), //subject olp q length
   st  char(4) default '', //olp type
   sct char(4) default '', //olp codon type
   sl  int(10), //length of homology olp
   evalue double,
   index qc(qc),
   index sc(sc),
   index qt(qt),
   index qcs(qc, qpg, qqg)) TYPE=MyISAM;


CREATE TABLE overlappingpairs (
   contig varchar(30) not null, 
   pg varchar(30) not null, 
   pg_product varchar(255) not null default '',
   pgs varchar(5), 
   qg varchar(30) not null,
   qg_product varchar(255) not null default '',
   qgs varchar(5), 
   type varchar(5) not null,
   codonType varchar(5),
   olpLength int(10),
   ps int(10),
   pe int(10),
   qs int(10),
   qe int(10),
   pgCDS text,
   polpAAs text,
   qgCDS text,
   qolpAAs text,
   index pg(pg),
   index qg(qg),
   index pqg(contig, pg, qg)) TYPE=MyISAM;

CREATE TABLE overlappingpairs (
   contig varchar(30) not null, //overlappen genes in the same contig
   pg varchar(30) not null, //p gene 
   pgs varchar(5), //p gene start
   qg varchar(30) not null, //q gene
   qgs varchar(5), //q gene start
   type varchar(5) not null,
   codonType varchar(5),
   olpLength int(10),
   ps int(10),
   pe int(10),
   qs int(10),
   qe int(10),
   pgCDS text,
   qgCDS text,
   index pg(pg),
   index qg(qg),
   index pqg(contig, pg, qg)) TYPE=MyISAM;

CREATE TABLE olp_domain(
   contig varchar(50) not null,
   pg varchar(10),
   pgs varchar(10),
   ps int(10),
   pe int(10),
   domain_pg tinytext,
   qg varchar(10),
   qgs varchar(5),
   qs int(10),
   qe int(10),
   domain_qg tinytext,
   olp_type varchar(5),
   olp_conType varchar(5),
   index pgc(contig, pg),
   index qgc(contig, qg)) TYPE=MyISAM;

CREATE TABLE mutualOlp(
   qc varchar(50) not null,
   qpg varchar(50),
   qqg varchar(50),
   sc varchar(50) not null,
   spg varchar(50),
   sqg varchar(50),
   qt varchar(6), 
   qct varchar(6),
   ql int(10),
   st varchar(6),
   sct varchar(6),
   sl int(10),
   index qc(qc),
   index sc(sc)) TYPE=MyISAM;


ALTER TAVLE olp_domain DROP INDEX accgene;
SHOW INDEX FROM olp_domain;

ALTER TABLE cdss MODIFY marker gpmarker int(5) AFTER gene_locus;
ALTER TABLE cdss ADD gqmarker int(5) AFTER gpmarker;
****************************************************************************************/

