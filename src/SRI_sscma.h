#ifdef SRI_SSCMA_EXPORTS
#define SRI_SSCMA_API __declspec(dllexport)
#else
#define SRI_SSCMA_API __declspec(dllimport)
#endif


extern "C" /* SRI_SSCMA_API */ 
int leaps(double *,int *,int *,int *,int *,int *,int *,int *,char **,
		  int *,int *,int *,double *,int *,int *,double *,double *,int *);


