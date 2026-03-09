#pragma once

#ifdef _MSC_VER
#define _CRT_SECURE_NO_WARNINGS
#endif


#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stddef.h>
#include <math.h>
#include <ctype.h>

#include "ModelicaUtilities.h"


typedef struct
{
	int order; // order of expansion
	int nrow; //number of rows
	int ncol; // number of columns
	int nq; //dimension of q
	int nqn;//dimension of higher order expansion
	int structure;// 0 =zero matrix, 1=diagonal, 2=symmetric, 3=no particular, 4 = unit matrix
	double* M0;// zero order matrix, nrow X ncol
	double* M1;// first order matrix, nrow X nq X ncol
	double* Mn;// higher order matrix
} taylor;

typedef struct
{
	int id;
	char rFrame[8];
	taylor orig;
	taylor phi;
	taylor psi;
	taylor AP;
} node;

typedef struct {
	int numNodes;
	int numModes;
	int useGeoStiffness;
	double mass;
	node* nodes;
	taylor mdCM;
	taylor J;
	taylor Ct;
	taylor Cr;
	taylor Me;
	taylor Gr;
	taylor Ge;
	taylor Oe;
	taylor ksigma;
	taylor Ke;
	taylor De;
} SID_Data;

typedef struct{
  SID_Data* sid;
} extObj;

// Index helper functions
//=================================

static inline int matrix2Index(int r, int c, int numR, int numC)
{
	return (r - 1) * numC + c - 1;
}

static inline int matrix3Index(int r, int q, int c, int numR, int numQ, int numC)
{
	return (q - 1) * numC * numR + (r - 1)*numC + c - 1;
}

//printing
//=================================

void printM_2(double* M, int dimRow, int dimCol) {
	int r =1;
	int c = 1;
	for (r = 1; r <= dimRow; r++) {
		for (c = 1; c <= dimCol; c++) {
			int i = matrix2Index(r, c, dimRow, dimCol);
			ModelicaFormatMessage("M[%d,%d]:%d = %f\n", r, c, i, M[i]);
		}
		ModelicaFormatMessage("\n");
	}
}

void printM_3(double* M, int dimRow, int dimQ, int dimCol) {
	int r = 1;
	int c = 1;
	int q = 1;
	for (r = 1; r <= dimRow; r++) {
		for (c = 1; c <= dimCol; c++) {
			for (q = 1; q <= dimQ; q++) {
				int i = matrix3Index(r, q, c, dimRow, dimQ, dimCol);
				ModelicaFormatMessage("M[%d,%d,%d]:%d = %f\n", r, q, c, i, M[i]);
			}
		}
	}
}

void printTaylor(taylor t) {
	ModelicaFormatMessage("\t order=%d\n", t.order);
	ModelicaFormatMessage("\tnrow=%d  ncol=%d  nq=%d  nqn=%d structure=%d \n", t.nrow, t.ncol, t.nq, t.nqn, t.structure);
	ModelicaFormatMessage("\t______________M0__________________\n");
	printM_2(t.M0, t.nrow, t.ncol);
	ModelicaFormatMessage("\t______________M1__________________\n");
	printM_3(t.M1, t.nrow, t.nq, t.ncol);
}


int getNextInteger(char* line, int startToken, int* value) {
	int len = strlen(line);
	int t = 0;
	int i = startToken;
	int found = 0;
	char intTokens[6] = "";
	char c;
	int v;
	while (i < len) {
		c = line[i];
		while (isdigit(c)) {
			strncat(intTokens, &c, 1);
			t++;
			i++;
			found = 1;
			c = line[i];
		}
		if (found) {
			break;
		}
		i++;
	}
	v = atoi(intTokens);
	*value = v;
	if (found)
		return i;
	else
		return -1;
}

int getNextReal(char* line, int startToken, double* value) {
	int len = strlen(line);
	int i = startToken;
	int found = 0;
	char baseToken[20] = "";
	char exponentToken[20] = "";
	char c;
	int seq = 0;
	double base;
	double exponent;
	double d;
	while (i < len) {
		c = line[i];
		//start real retrieval
		if (isdigit(c) || c == 45) {
			seq = 1;
		}
		while ((0 < seq) && (seq< 3)) {
			//some digit
			if (isdigit(c) || c == 45 || c == 46) {
				//post comma
				if (seq == 1) {
					strncat(baseToken, &c, 1);
				}
				else if (seq == 2) {
					strncat(exponentToken, &c, 1);
				}
			}
			//exponent
			if (c == 'D' && seq == 1) {
				seq = 2;

			}
			if (c == '\n') {
				break;
			}
			i++;
			found = 1;
			c = line[i];
		}
		i++;
	}
	base = atof(baseToken);
	exponent = atof(exponentToken);
	d = base * pow(10, exponent);
	*value = d;
	if (found)
		return i;
	else
		return -1;
}


void parseFirstLine(char* line, SID_Data* sid) {
	int startToken = 0;
	int* numNodes = calloc(sizeof(int), 1);
	int* numModes = calloc(sizeof(int), 1);
	startToken = getNextInteger(line, startToken, numNodes);
	startToken = getNextInteger(line, startToken, numModes);
	sid->numNodes = numNodes[0];
	sid->numModes = numModes[0];
	free(numNodes);
	free(numModes);
}

void parseRefMod(char* line, int bufferSize, FILE* file, SID_Data* sid) {
	while (fgets(line, bufferSize, file) != NULL) {
		if (strstr(line, " mass ")) {
			int startToken = 0;
			double mass = 0;
			int found = getNextReal(line, startToken, &mass);
			sid->mass = mass;
		}
		else if (strstr(line, "end refmod"))
		{
			return;
		}
	}
	ModelicaFormatMessage("Did not found end of refmod\n");
}

taylor initTaylor(int order, int nr, int ncol, int nq, int nqn, int structure) {
	taylor t;// = { .order = 0,.nrow = nr,.ncol = ncol,.nq = nq,.nqn = nqn,.structure = structure,.M0 = NULL,.M1 = NULL,.Mn = NULL };
	t.order = 0;
	t.nrow = nr;
	t.ncol = ncol;
	t.nq = nq;
	t.nqn = nqn;
	t.structure = structure;
	t.Mn = NULL;
	t.M0 = calloc(t.nrow*t.ncol, sizeof(double));
	t.M1 = calloc(t.nrow*t.nq*t.ncol, sizeof(double));
	if (t.nqn != 0)
		t.Mn = calloc(t.nrow*t.nqn*t.ncol, sizeof(double));
	return t;
}

taylor parseTaylor(char* line, int bufferSize, FILE* file) {
	int order = 0;
	int nrow = 0;
	int ncol = 0;
	int nq = 0;
	int nqn = 0;
	int structure = 0;
	taylor t;
	int gotObjects = 0;
	while (fgets(line, bufferSize, file) != NULL) {
		if (gotObjects == 6) {
			t = initTaylor(order, nrow, ncol, nq, nqn, structure);
			gotObjects = 7;
		}

		if (strstr(line, "order ")) {
			getNextInteger(line, 0, &order);
			gotObjects++;
		}
		else if (strstr(line, "nrow ")) {
			getNextInteger(line, 0, &nrow);
			gotObjects++;
		}
		else if (strstr(line, "ncol ")) {
			getNextInteger(line, 0, &ncol);
			gotObjects++;
		}
		else if (strstr(line, "nqn ")) {
			getNextInteger(line, 0, &nqn);
			gotObjects++;
		}
		else if (strstr(line, "nq ")) {
			getNextInteger(line, 0, &nq);
			gotObjects++;
		}
		else if (strstr(line, "structure ")) {
			getNextInteger(line, 0, &structure);
			gotObjects++;
		}
		else if (strstr(line, "m0")) {
			int idx1 = 0;
			int idx2 = 0;
			int tokenIdx = 0;
			double val = 0.0;
			int i;
			tokenIdx = getNextInteger(line, tokenIdx, &idx1);// this is the 0 in m0
			tokenIdx = getNextInteger(line, tokenIdx + 1, &idx1);
			tokenIdx = getNextInteger(line, tokenIdx + 1, &idx2);
			getNextReal(line, tokenIdx, &val);
			i = matrix2Index(idx1, idx2, nrow, ncol);
			//ModelicaFormatMessage("set M0 i %d  = %f\n",i,val);
			//ModelicaFormatMessage("set M0 i %d  = %f\n",i,val);
			t.M0[i] = val;
			//printM_2(t.M0, nrow, ncol);
		}
		else if (strstr(line, "m1")) {
			int idx1 = 0;
			int idx2 = 0;
			int idx3 = 0;
			int i = 0;
			int tokenIdx = 0;
			double val = 0.0;
			tokenIdx = getNextInteger(line, tokenIdx, &idx1);// this is the 0 in m0
			tokenIdx = getNextInteger(line, tokenIdx + 1, &idx1);
			tokenIdx = getNextInteger(line, tokenIdx + 1, &idx2);
			tokenIdx = getNextInteger(line, tokenIdx + 1, &idx3);
			getNextReal(line, tokenIdx, &val);
			i = matrix3Index(idx1, idx2, idx3, nrow, nq, ncol);
			//ModelicaFormatMessage("set M1 i %d  = %f\n", i, val);
			t.M1[i] = val;
		}
		else if (strstr(line, "end "))
		{
			//printTaylor(t);
			return t;
		}
	}
	ModelicaFormatMessage("Did not found end of taylor\n");
	return t;
}

void parseNode(char* line, int bufferSize, FILE* file, SID_Data* sid, int nodeIdx) {
	taylor t;
	//ModelicaFormatMessage("ParseNode  %d \n",nodeIdx);
	while (fgets(line, bufferSize, file) != NULL) {
		if (strstr(line, "origin")) {
			//ModelicaFormatMessage("origin\n");
			t = parseTaylor(line, bufferSize, file);
			sid->nodes[nodeIdx].orig = t;
		}
		if (strstr(line, "phi")) {
			//ModelicaFormatMessage("phi\n");
			t = parseTaylor(line, bufferSize, file);
			sid->nodes[nodeIdx].phi = t;
		}
		if (strstr(line, "psi")) {
			//ModelicaFormatMessage("psi\n");
			t = parseTaylor(line, bufferSize, file);
			sid->nodes[nodeIdx].psi = t;
		}
		if (strstr(line, "AP")) {
			//ModelicaFormatMessage("AP\n");
			t = parseTaylor(line, bufferSize, file);
			sid->nodes[nodeIdx].AP = t;
		}
		else if (strstr(line, "end node"))
		{
			return;
		}
	}
	ModelicaFormatMessage("Did not found end of node\n");
}

void parseNodes(char* line, int bufferSize, FILE* file, SID_Data* sid) {
	int nodeIdx = 0;
	while (fgets(line, bufferSize, file) != NULL) {
		if (strstr(line, "new node")) {
			int* nodeID = calloc(sizeof(int), 1);
			getNextInteger(line, 0, nodeID);
			sid->nodes[nodeIdx].id = nodeID[0];
			parseNode(line, bufferSize, file, sid, nodeIdx);
			nodeIdx++;
		}
		else if (strstr(line, "end node"))
		{
			if (nodeIdx == sid->numNodes) {
			}
		}
		else if (strstr(line, "end frame")) {
			return;
		}
	}
	ModelicaFormatMessage("Did not found end of frame\n");
}

taylor getNodeTaylorByName(node* node, const char* name)
{
	//ModelicaFormatMessage("getNodeTaylorByName %s\n",name);
	if(!strcmp(name,"origin")){
		return node->orig;
	}
	else if (!strcmp(name,"phi")){
	    return node->phi;
	}
	else if (!strcmp(name,"psi")){
		return node->psi;
	}
	else if (!strcmp(name,"AP")){
		return node->AP;
	}
	else{
	ModelicaFormatMessage("getNodeTaylorByName: the given name is not defined: %s\n",name);
	return node->orig;
	}
}


void freeTaylor(taylor* t){
	free(t->M0);
	free(t->M1);
	free(t->Mn);
	}


void SIDFileDestructor_C(void* p_sid) {
	/*
	SID_Data* sid = (SID_Data*)p_sid;
	for (int n=0;n<sid->numNodes;n++){
		freeTaylor(&sid->nodes[n].orig);
		freeTaylor(&sid->nodes[n].phi);
		freeTaylor(&sid->nodes[n].psi);
		freeTaylor(&sid->nodes[n].AP);
	}
	free(sid->nodes);
	freeTaylor(&sid->mdCM);
	freeTaylor(&sid->J);
	freeTaylor(&sid->Ct);
	freeTaylor(&sid->Cr);
	freeTaylor(&sid->Me);
	freeTaylor(&sid->Gr);
	freeTaylor(&sid->Ge);
	freeTaylor(&sid->Oe);
	freeTaylor(&sid->ksigma);
	freeTaylor(&sid->Ke);
	freeTaylor(&sid->De);
	free(sid);
	*/
}

void* SIDFileConstructor_C(const char* fileName)
{
	FILE *sidFile;
	SID_Data* sid;
    int bufferSize;
	char *line;
	int init;
   	taylor t;

	if ((sidFile = fopen(fileName, "r")) == NULL) {
		ModelicaFormatMessage("Could not open %s\n", fileName);
	}
	sid = calloc(1, sizeof(SID_Data));
	ModelicaFormatMessage("Allocated SID struct of size %zu\n", sizeof(SID_Data));

	//read line by line ( http://openbook.rheinwerk-verlag.de/c_von_a_bis_z/016_c_ein_ausgabe_funktionen_016.htm#mjcea47bd6d32a4a8f51be329a672845d7 )
	bufferSize = 512;
	line = malloc(sizeof(char) * 512);

	init = 0;
	while (fgets(line, bufferSize, sidFile) != NULL) {
		if (!init) {
			parseFirstLine(line, sid);
			init = 1;
			sid->nodes = calloc(sizeof(node), sid->numNodes);
			//skip the rest of the first line
			fgets(line, bufferSize, sidFile); 
		}
		if (strstr(line, "refmod")) {
			parseRefMod(line, bufferSize, sidFile, sid);
		}
		else if (strstr(line, "frame")) {
			parseNodes(line, bufferSize, sidFile, sid);
		}
		else if (strstr(line, "mdCM")) {
			t = parseTaylor(line, bufferSize, sidFile);
			sid->mdCM = t;
		}
		else if (strstr(line, "J")) {
			t = parseTaylor(line, bufferSize, sidFile);
			sid->J = t;
		}
		else if (strstr(line, "Ct")) {
			t = parseTaylor(line, bufferSize, sidFile);
			sid->Ct = t;
		}
		else if (strstr(line, "Cr")) {
			t = parseTaylor(line, bufferSize, sidFile);
			sid->Cr = t;
		}
		else if (strstr(line, "Me")) {
			t = parseTaylor(line, bufferSize, sidFile);
			sid->Me = t;
		}
		else if (strstr(line, "Gr")) {
			t = parseTaylor(line, bufferSize, sidFile);
			sid->Gr = t;
		}
		else if (strstr(line, "Ge")) {
			t = parseTaylor(line, bufferSize, sidFile);
			sid->Ge = t;
		}
		else if (strstr(line, "Oe")) {
			t = parseTaylor(line, bufferSize, sidFile);
			sid->Oe = t;
		}
		else if (strstr(line, "ksigma")) {
			t = parseTaylor(line, bufferSize, sidFile);
			sid->ksigma = t;
		}
		else if (strstr(line, "Ke")) {
			t = parseTaylor(line, bufferSize, sidFile);
			sid->Ke = t;
		}
		else if (strstr(line, "De")) {
			t = parseTaylor(line, bufferSize, sidFile);
			sid->De = t;
		}
	}
	fclose(sidFile);
	ModelicaFormatMessage("Read SID FIle %s , mass is %f\n",fileName,sid->mass);
	return (void*)sid;
}


//retrieve basic data
//=================================

double getMass(void* p_sid){
	SID_Data* sid = (SID_Data*)p_sid;
	double	d = sid->mass;
	//ModelicaFormatMessage("getMass %f     %f \n",sid->mass,d);
	return d;
}



//Taylor Retrieval
//=================================

taylor getTaylorByName(SID_Data* sid, const char* name)
{
	//ModelicaFormatMessage("getTaylorByName %s\n",name);
	if(!strcmp(name,"mdCM")){
		return sid->mdCM;
	}
	else if (!strcmp(name,"J")){
	    return sid->J;
	}
	else if (!strcmp(name,"Ct")){
		return sid->Ct;
	}
	else if (!strcmp(name,"Cr")){
		return sid->Cr;
	}
	else if (!strcmp(name,"Me")){
		return sid->Me;
	}
	else if (!strcmp(name,"Gr")){
		return sid->Gr;	
	}
	else if (!strcmp(name,"Ge")){
		return sid->Ge;
	}
	else if (!strcmp(name,"Oe")){
		return sid->Oe;
	}
	else if (!strcmp(name,"ksigma")){
		return sid->ksigma;
	}
	else if (!strcmp(name,"Ke")){
		return sid->Ke;
	}
	else if (!strcmp(name,"De")){
		return sid->De;
	}
	else{
	ModelicaFormatMessage("getTaylorByName: the given name is not defined: %s\n",name);
	return sid->mdCM;
	}
}

//M0 for modal objects
void getM0(void* p_sid, const char* taylorName, double* m0, size_t nr, size_t nc){
	SID_Data* sid = (SID_Data*)p_sid;
	taylor t = getTaylorByName((sid), taylorName);
	if(!((int)nr==t.nrow && (int)nc==t.ncol)){
		ModelicaFormatMessage(" getM0: %s the given dimensions [%zu, %zu] are not equal to the stored matrix dimension [%d %d]\n",taylorName, nr, nc, t.nrow,t.ncol);
	}
	else{
		int s = sizeof(double)*(int)nr*(int)nc;
	    memcpy(m0,t.M0,s);
	}
}
//M1 for modal objects
void getM1(void* p_sid, const char* taylorName, double* m1, size_t nr, size_t nq, size_t nc){
	SID_Data* sid = (SID_Data*)p_sid;
	//ModelicaFormatMessage("getM1 %s [%d, %d,  %d]\n",taylorName,nr,nq,nc);
	taylor t = getTaylorByName((sid), taylorName);

	if((int)nr!=t.nrow || (int)nc!=t.ncol || (int)nq!=t.nq){
		ModelicaFormatMessage(" getM1: the given dimensions [%zu, %zu, %zu] are not equal to the stored matrix dimension [%d, %d, %d]\n",nr,nq,nc, t.nrow,t.nq,t.ncol);
	}
	else{
	  int s = sizeof(double)*(int)nr*(int)nc*(int)nq;
	  memcpy(m1,t.M1,s);
	}
}

//M0 for node object
void getM0Node(void* p_sid, const char* taylorName, int nodeIdx, double* m0, size_t nr, size_t nc){
	SID_Data* sid = (SID_Data*)p_sid;
	node n;
	taylor t;
	if (nodeIdx > sid->numNodes){
		ModelicaFormatMessage(" getM0Node: the given nodeidx  [%d] is bigger than the number of nodes [%d]\n", nodeIdx, sid->numNodes);
	}
	n = sid->nodes[nodeIdx-1];
	t = getNodeTaylorByName(&n, taylorName);
	if((int)nr!=t.nrow || (int)nc!=t.ncol){
		ModelicaFormatMessage(" getM0Node: the given dimensions [%zu, %zu] are not equal to the stored matrix dimension [%d %d]\n", nr, nc, t.nrow,t.ncol);
	}
	else{
	  int s = sizeof(double)*(int)nr*(int)nc;
	  memcpy(m0,t.M0,s);
	}
	
	
}
//M1 for node object
void getM1Node(void* p_sid, const char* taylorName, int nodeIdx, double* m1, size_t nr, size_t nq, size_t nc){
	SID_Data* sid = (SID_Data*)p_sid;
	node n;
	taylor t;
	if (nodeIdx > sid->numNodes){
		ModelicaFormatMessage(" getM0Node: the given nodeidx  [%d] is bigger than the number of nodes [%d]\n", nodeIdx, sid->numNodes);
	}
	n = sid->nodes[nodeIdx-1];
	t = getNodeTaylorByName(&n, taylorName);
			
	if((int)nr!=t.nrow || (int)nc!=t.ncol || (int)nq!=t.nq){
		ModelicaFormatMessage(" getM1: the given dimensions [%zu, %zu, %zu] are not equal to the stored matrix dimension [%d, %d, %d]\n",nr,nq,nc, t.nrow,t.nq,t.ncol);
	}
	else{
	  int s = sizeof(double)*(int)nr*(int)nc*(int)nq;
	  memcpy(m1,t.M1,s);
	}

}
