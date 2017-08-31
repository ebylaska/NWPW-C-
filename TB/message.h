#ifndef _MESSAGE_H_
#define _MESSAGE_H_


extern "C" {
void Header_Message();
void Header_TBMD_Message();
void Input_Data_Message(const int,    const int,
			       const double, const double,
			       const int, char**, double**,
			       const char*);
void Input_Data_TBMD_Message(const int,    const int,
			       const double, const double,
			       const int, char**, double**, double**,
			       const double, const double,
			       const double, const double,
			       const char*);
void Output_Data_Message(const double, 
                                const double, 
                                const double,
                                const double,
                                const int, char**, double**);
void Output_Data_TBMD_Message(const double, 
                                const double, 
                                const double,
                                const double,
                                const int, char**, double**, double**);

void Start_Iteration_Message();
void Iteration_Message(const int, const int, const double*);
void End_Iteration_Message(const int);
void Timing_Message(const double, const double,
                           const double, const double,
			   const double);
void Start_Eigenvalue_Message(const int);
void Eigenvalue_Message(const int, const int, const double);

}
#endif
