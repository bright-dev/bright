// Header for general library file.

#if !defined(_bright_)
#define _bright_

//standard libraries
#include <string>
#include <string.h>
#include <sstream>
#include <cctype>
#include <stdlib.h>
#include <iostream>
#include <math.h>

/*** Macros ***/
#define length_array(a) ( sizeof ( a ) / sizeof ( *a ) )

namespace bright {

//Bright Globals
void BrightStart ();

extern std::string BRIGHT_DATA;

//String Transformations
static std::string alphabet = "ABCDEFGHIJKLMNOPQRSTUVWXYZ";
static std::string digits = "0123456789";

std::string to_str (int);
std::string to_str (double);
std::string to_str (bool);

int to_int (std::string);

double to_dbl (std::string);

std::string ToUpper(std::string);

std::string ToLower(std::string);

std::string getFlag(char [], int);

std::string Strip(std::string, std::string);

std::string MultiStrip(std::string, std::string);

std::string LastChar(std::string);

std::string SubFromEnd(std::string, int = -1, int = 1);

bool ChainGreaterCompare(int, int, int);

bool SubInString(std::string, std::string);

//Array Methods
template <class T>
int find_index(T val, T * arr, int arr_len = -1)
{
    //Finds an element 'val' in array 'arr'
    //returns the index of val's first location
    //returns -1 if not found.

    if (arr_len < 0)
        arr_len = sizeof(arr) / sizeof(T);

    for (int n = 0; n < arr_len; n++)
    {
        if (val == arr[n])
            return n;
    };

    return -1;
};

int find_index_char(char *, char **, int = -1);

//Math Helpers
extern const double pi;
extern const double N_A;	//Avagardo's Number
extern const double bpcm2; 	//barns per cm^2

double slope (double, double, double, double);
double SolveLine (double, double, double, double, double);

double TANH(double);
double COTH(double);

};

#endif
