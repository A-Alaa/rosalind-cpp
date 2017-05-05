#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* as PatternCount.c works this function instead return an array of starting positions where the subtext occurs in the main text and considering working with large texts file I/O stream used instead of passing the large text to "char** argv" */


int main(int argc , char** argv ){
	
	FILE* fp;
	fp=fopen(argv[1],"r");
	
	int kmer=strlen(argv[2]);
	char* temp=(char*)malloc(kmer);
	
	int posArraySize=0;
	int* posArray=(int*)realloc(NULL,sizeof(int)*posArraySize);

	int i=0;
	while(1){
		fseek(fp,i,SEEK_SET);
		if(fread(temp,sizeof(char),kmer,fp)<kmer) break;
		if(strncmp(temp,argv[2],kmer)==0){
			posArray=(int*)realloc(posArray,sizeof(int)*(posArraySize+1));
			posArray[posArraySize++]=i;
		}
		i++;
	}
	
	for(i=0;i<posArraySize;i++) printf("%d ",posArray[i]);
	printf("\n");
	
}		
