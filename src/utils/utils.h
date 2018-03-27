#ifndef __INFOLI_UTILS_H__
#define __INFOLI_UTILS_H__

typedef float mod_prec;

void print_usage();
void read_parameters(char *paramsFileName, mod_prec *values);
void removeSubstring(char *s, const char *toremove);
void stopAtSubstring(char *s, const char *toremove);
int ReadFileLine(FILE *pInFile, mod_prec *iAppArray);

static inline mod_prec min(mod_prec a, mod_prec b){

	return (a < b) ? a : b;
}


#endif /* __INFOLI_UTILS_H__ */
