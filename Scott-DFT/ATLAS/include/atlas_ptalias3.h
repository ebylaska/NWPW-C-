#ifndef ATLAS_PTALIAS3_H
#define ATLAS_PTALIAS3_H
/*
 * Real BLAS
 */
   #ifdef ATL_ANTOINE_THREADS
      #define ATL_ssyr2k  ATL_sptsyr2k
   #elif defined(ATL_OMP_THREADS)
      #define ATL_ssyr2k  ATL_stompsyr2k
   #else
      #define ATL_ssyr2k  ATL_stsyr2k
   #endif
   #ifdef ATL_ANTOINE_THREADS
      #define ATL_ssyrk   ATL_sptsyrk
   #elif defined(ATL_OMP_THREADS)
      #define ATL_ssyrk   ATL_stompsyrk
   #else
      #define ATL_ssyrk   ATL_stsyrk
   #endif
   #ifdef ATL_ANTOINE_THREADS
      #define ATL_ssymm   ATL_sptsymm
   #elif defined(ATL_OMP_THREADS)
      #define ATL_ssymm   ATL_stompsymm
   #else
      #define ATL_ssymm   ATL_stsymm
   #endif
   #ifdef ATL_ANTOINE_THREADS
      #define ATL_strmm   ATL_spttrmm
   #elif defined(ATL_OMP_THREADS)
      #define ATL_strmm   ATL_stomptrmm
   #else
      #define ATL_strmm   ATL_sttrmm
   #endif
   #ifdef ATL_ANTOINE_THREADS
      #define ATL_strsm   ATL_spttrsm
   #elif defined(ATL_OMP_THREADS)
      #define ATL_strsm   ATL_stomptrsm
   #else
      #define ATL_strsm   ATL_sttrsm
   #endif
   #ifdef ATL_ANTOINE_THREADS
      #define ATL_sgemm   ATL_sptgemm
   #elif defined(ATL_OMP_THREADS)
      #define ATL_sgemm   ATL_stompgemm
   #else
      #define ATL_sgemm   ATL_stgemm
   #endif

   #ifdef ATL_ANTOINE_THREADS
      #define ATL_dsyr2k  ATL_dptsyr2k
   #elif defined(ATL_OMP_THREADS)
      #define ATL_dsyr2k  ATL_dtompsyr2k
   #else
      #define ATL_dsyr2k  ATL_dtsyr2k
   #endif
   #ifdef ATL_ANTOINE_THREADS
      #define ATL_dsyrk   ATL_dptsyrk
   #elif defined(ATL_OMP_THREADS)
      #define ATL_dsyrk   ATL_dtompsyrk
   #else
      #define ATL_dsyrk   ATL_dtsyrk
   #endif
   #ifdef ATL_ANTOINE_THREADS
      #define ATL_dsymm   ATL_dptsymm
   #elif defined(ATL_OMP_THREADS)
      #define ATL_dsymm   ATL_dtompsymm
   #else
      #define ATL_dsymm   ATL_dtsymm
   #endif
   #ifdef ATL_ANTOINE_THREADS
      #define ATL_dtrmm   ATL_dpttrmm
   #elif defined(ATL_OMP_THREADS)
      #define ATL_dtrmm   ATL_dtomptrmm
   #else
      #define ATL_dtrmm   ATL_dttrmm
   #endif
   #ifdef ATL_ANTOINE_THREADS
      #define ATL_dtrsm   ATL_dpttrsm
   #elif defined(ATL_OMP_THREADS)
      #define ATL_dtrsm   ATL_dtomptrsm
   #else
      #define ATL_dtrsm   ATL_dttrsm
   #endif
   #ifdef ATL_ANTOINE_THREADS
      #define ATL_dgemm   ATL_dptgemm
   #elif defined(ATL_OMP_THREADS)
      #define ATL_dgemm   ATL_dtompgemm
   #else
      #define ATL_dgemm   ATL_dtgemm
   #endif

/*
 * Complex BLAS
 */
   #ifdef ATL_ANTOINE_THREADS
      #define ATL_ctrmm     ATL_cpttrmm
   #elif defined(ATL_OMP_THREADS)
      #define ATL_ctrmm     ATL_ctomptrmm
   #else
      #define ATL_ctrmm     ATL_cttrmm
   #endif
   #ifdef ATL_ANTOINE_THREADS
      #define ATL_cher2k    ATL_cpther2k
   #elif defined(ATL_OMP_THREADS)
      #define ATL_cher2k    ATL_ctompher2k
   #else
      #define ATL_cher2k    ATL_cther2k
   #endif
   #ifdef ATL_ANTOINE_THREADS
      #define ATL_csyr2k    ATL_cptsyr2k
   #elif defined(ATL_OMP_THREADS)
      #define ATL_csyr2k    ATL_ctompsyr2k
   #else
      #define ATL_csyr2k    ATL_ctsyr2k
   #endif
   #ifdef ATL_ANTOINE_THREADS
      #define ATL_cherk     ATL_cptherk
   #elif defined(ATL_OMP_THREADS)
      #define ATL_cherk     ATL_ctompherk
   #else
      #define ATL_cherk     ATL_ctherk
   #endif
   #ifdef ATL_ANTOINE_THREADS
      #define ATL_csyrk     ATL_cptsyrk
   #elif defined(ATL_OMP_THREADS)
      #define ATL_csyrk     ATL_ctompsyrk
   #else
      #define ATL_csyrk     ATL_ctsyrk
   #endif
   #ifdef ATL_ANTOINE_THREADS
      #define ATL_chemm     ATL_cpthemm
   #elif defined(ATL_OMP_THREADS)
      #define ATL_chemm     ATL_ctomphemm
   #else
      #define ATL_chemm     ATL_cthemm
   #endif
   #ifdef ATL_ANTOINE_THREADS
      #define ATL_csymm     ATL_cptsymm
   #elif defined(ATL_OMP_THREADS)
      #define ATL_csymm     ATL_ctompsymm
   #else
      #define ATL_csymm     ATL_ctsymm
   #endif
   #ifdef ATL_ANTOINE_THREADS
      #define ATL_cgemm     ATL_cptgemm
   #elif defined(ATL_OMP_THREADS)
      #define ATL_cgemm     ATL_ctompgemm
   #else
      #define ATL_cgemm     ATL_ctgemm
   #endif
   #ifdef ATL_ANTOINE_THREADS
      #define ATL_ctrsm     ATL_cpttrsm
   #elif defined(ATL_OMP_THREADS)
      #define ATL_ctrsm     ATL_ctomptrsm
   #else
      #define ATL_ctrsm     ATL_cttrsm
   #endif

   #ifdef ATL_ANTOINE_THREADS
      #define ATL_ztrmm     ATL_zpttrmm
   #elif defined(ATL_OMP_THREADS)
      #define ATL_ztrmm     ATL_ztomptrmm
   #else
      #define ATL_ztrmm     ATL_zttrmm
   #endif
   #ifdef ATL_ANTOINE_THREADS
      #define ATL_zher2k    ATL_zpther2k
   #elif defined(ATL_OMP_THREADS)
      #define ATL_zher2k    ATL_ztompher2k
   #else
      #define ATL_zher2k    ATL_zther2k
   #endif
   #ifdef ATL_ANTOINE_THREADS
      #define ATL_zsyr2k    ATL_zptsyr2k
   #elif defined(ATL_OMP_THREADS)
      #define ATL_zsyr2k    ATL_ztompsyr2k
   #else
      #define ATL_zsyr2k    ATL_ztsyr2k
   #endif
   #ifdef ATL_ANTOINE_THREADS
      #define ATL_zherk     ATL_zptherk
   #elif defined(ATL_OMP_THREADS)
      #define ATL_zherk     ATL_ztompherk
   #else
      #define ATL_zherk     ATL_ztherk
   #endif
   #ifdef ATL_ANTOINE_THREADS
      #define ATL_zsyrk     ATL_zptsyrk
   #elif defined(ATL_OMP_THREADS)
      #define ATL_zsyrk     ATL_ztompsyrk
   #else
      #define ATL_zsyrk     ATL_ztsyrk
   #endif
   #ifdef ATL_ANTOINE_THREADS
      #define ATL_zhemm     ATL_zpthemm
   #elif defined(ATL_OMP_THREADS)
      #define ATL_zhemm     ATL_ztomphemm
   #else
      #define ATL_zhemm     ATL_zthemm
   #endif
   #ifdef ATL_ANTOINE_THREADS
      #define ATL_zsymm     ATL_zptsymm
   #elif defined(ATL_OMP_THREADS)
      #define ATL_zsymm     ATL_ztompsymm
   #else
      #define ATL_zsymm     ATL_ztsymm
   #endif
   #ifdef ATL_ANTOINE_THREADS
      #define ATL_zgemm     ATL_zptgemm
   #elif defined(ATL_OMP_THREADS)
      #define ATL_zgemm     ATL_ztompgemm
   #else
      #define ATL_zgemm     ATL_ztgemm
   #endif
   #ifdef ATL_ANTOINE_THREADS
      #define ATL_ztrsm     ATL_zpttrsm
   #elif defined(ATL_OMP_THREADS)
      #define ATL_ztrsm     ATL_ztomptrsm
   #else
      #define ATL_ztrsm     ATL_zttrsm
   #endif

#endif
