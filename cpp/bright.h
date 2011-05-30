// Header for general library file.

#if !defined(_bright_)
#define _bright_

//standard libraries
#include <string>
#include <string.h>
#include <sstream>
#include <iostream>
#include <cctype>
#include <stdlib.h>
#include <iostream>
#include <math.h>
#include <exception>
#include <sys/stat.h> 
#include <vector>
#include <algorithm>
#include <typeinfo>

/*** Macros ***/
#define length_array(a) ( sizeof ( a ) / sizeof ( *a ) )

#ifdef _WIN32
    #define isnan(x) ((x) != (x))
#endif

namespace bright {

//Bright Globals
void bright_start ();

extern std::string BRIGHT_DATA;

//String Transformations
static std::string alphabet = "ABCDEFGHIJKLMNOPQRSTUVWXYZ";
static std::string digits = "0123456789";
static std::string words = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789_";


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

std::string StrictReplaceAll(std::string, std::string, std::string);

std::string LastChar(std::string);

std::string SubFromEnd(std::string, int = -1, int = 1);

bool ChainGreaterCompare(int, int, int);

bool SubInString(std::string, std::string);

std::string natural_naming(std::string);


// Vectorized Functions
std::vector<double> delta_vector(double, std::vector<double>);
std::vector<double> normalized_delta(double, std::vector<double>);

bool sorted_index_comparator(std::pair<int, double>, std::pair<int, double>);
std::vector<int> sorted_index(std::vector<double>);

std::vector<double> y_x_factor_interpolation(double, std::vector<double>, std::vector<double>);

std::vector< std::vector<double> > vector_outer_product(std::vector<double>, std::vector<double>);
std::vector< std::vector<double> > matrix_inverse(std::vector< std::vector<double> >);
std::vector< std::vector<double> > matrix_addition(std::vector< std::vector<double> >, std::vector< std::vector<double> >);
std::vector< std::vector<double> > matrix_multiplication(std::vector< std::vector<double> >, std::vector< std::vector<double> >);

std::vector< std::vector<double> > scalar_matrix_product(double, std::vector< std::vector<double> >);
std::vector<double> scalar_matrix_vector_product(double, std::vector< std::vector<double> >, std::vector<double>);

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
extern const double barns_per_cm2; 	// barns per cm^2
extern const double cm2_per_barn; 	// cm^2 per barn
extern const double sec_per_day; 	// seconds per day

double slope (double, double, double, double);
double SolveLine (double, double, double, double, double);

double TANH(double);
double COTH(double);


// File Helpers
bool FileExists(std::string); 



/******************/
/*** Exceptions ***/
/******************/


class FileNotFound : public std::exception
{
public:
    FileNotFound () {};

    ~FileNotFound () throw () {};

    FileNotFound(std::string fname)
    {
            filename = fname;
    };

    virtual const char* what() const throw()
    {
        std::string FNFstr ("File not found: ");
        if (!filename.empty())
            FNFstr += filename;

        return (const char *) FNFstr.c_str();
    };

private:
    std::string filename;
};




class BadFuelForm : public std::exception
{
//Exception for valid fuel form.
public:
    BadFuelForm () {};
    ~BadFuelForm () throw () {};

    static char * name ()
    {
        return (char *) "BadFuelForm";
    };

    virtual const char* what() const throw()
    {
        std::string BFFstr ("FUEL COMPOSITION NOT COMPUTABLE!");
        return (const char *) BFFstr.c_str();
    };
};





class VectorSizeError : public std::exception
{
//Exception for valid fuel form.
public:
    VectorSizeError () {};
    ~VectorSizeError () throw () {};

    static char * name ()
    {
        return (char *) "VectorSizeError";
    };

    virtual const char* what() const throw()
    {
        std::string VWSstr ("Vector is of the wrong size.");
        return (const char *) VWSstr.c_str();
    };
};





class BisectionMethodNotPerformed : public std::exception
{
//Exception for when the bisection method is not calculated.
public:
    BisectionMethodNotPerformed ()
    {
        errstr = "Bisection method was not performed.";
    };
    BisectionMethodNotPerformed (std::string calctype)
    {
        errstr = "Bisection method durring " + calctype + " calculation was not performed.";
    };
    ~BisectionMethodNotPerformed () throw () {};

    static char * name ()
    {
        return (char *) "BisectionMethodNotPerformed";
    };

    virtual const char* what() const throw()
    {
        return (const char *) errstr.c_str();
    };
private:
    std::string errstr;
};




template <class T>
class sparse_matrix_entry
{
public:
    int row;
    int col;
    T val;

    sparse_matrix_entry() {};
    ~sparse_matrix_entry() {};

    sparse_matrix_entry(int i, int j, T v)
    {
        row = i;
        col = j;
        val = v;
    };
};


template <class T>
bool cmp_by_row (sparse_matrix_entry<T> a,sparse_matrix_entry<T> b) 
{
    if (a.row != b.row)
        return (a.row < b.row);
    else
        return (a.col < b.col);
};


template <class T>
bool cmp_by_col (sparse_matrix_entry<T> a,sparse_matrix_entry<T> b) 
{
    if (a.col != b.col)
        return (a.col < b.col);
    else
        return (a.row < b.row);
};


template<class InputIterator, class T>
InputIterator find_row( InputIterator first, InputIterator last, const T& value )
{
    for ( ;first!=last; first++) 
    {
        if ((*first).row == value) 
            break;

        if (value < (*first).row)
            first = last - 1;
    };

    return first;
};



template<class InputIterator, class T>
InputIterator find_col( InputIterator first, InputIterator last, const T& value )
{
    for ( ;first!=last; first++) 
    {
        if ((*first).col == value) 
            break;

        if (value < (*first).col) 
            first = last - 1 ;
    };

    return first;
};



template <class T>
class SparseMatrix
{
public:
    int nrows, ncols;
    std::vector< sparse_matrix_entry<T> > sm;

    SparseMatrix(){};
    ~SparseMatrix(){};

    SparseMatrix(int N, int nr=0, int nc = 0)
    {
        nrows = nr;
        ncols = nc;

        sm = std::vector< sparse_matrix_entry <T> >();
        sm.reserve(N);
    };


    int size()
    {
        return sm.size();
    };


    void push_back(int i, int j, T value)
    {
        sm.push_back(sparse_matrix_entry<T>(i, j, value));
    };


    T at(int i, int j)
    {
        typename std::vector< sparse_matrix_entry<T> >::iterator a_iter = find_row(sm.begin(), sm.end(), i);
        while (i == (*a_iter).row)
        {
            if (j == (*a_iter).col)
                return (*a_iter).val;

            a_iter++;
        };
        return 0.0;
    };


    void sort_by_row()
    {
        std::sort(sm.begin(), sm.end(), cmp_by_row<T>);
    };


    void sort_by_col()
    {
        std::sort(sm.begin(), sm.end(), cmp_by_col<T>);
    };


    void clean_up()
    {
        // First, get all of your ducks in a row
        sort_by_row();

        int n, N;
        N = sm.size();
        std::vector<int> bad_ind = std::vector<int>();

        // Calculate indices to remove
        for (n = N - 1; 0 < n; n--)
        {
            if ((sm[n].row == sm[n-1].row) && (sm[n].col == sm[n-1].col))
                bad_ind.push_back(n);
            else if (sm[n].val == 0.0)
                bad_ind.push_back(n);
        }; 

        // remove the offending indices
        int p, P;
        P = bad_ind.size();
        for (p = 0; p < P; p++)
            sm.erase(sm.begin()+bad_ind[p]);

        // Save some space
        sm.resize(sm.size());
    };


    friend std::ostream& operator<< (std::ostream& out, SparseMatrix<T> & A) 
    {
        int n = 0;
        int N = A.size();
        
        out << "Sparse Matrix [" << A.nrows << ", " << A.ncols << "] (" << N << ")\n";
        for (n = 0; n < N; n++)
            out << "  (" << A.sm[n].row << ", " << A.sm[n].col << ") = " << A.sm[n].val << "\n";

        return out;
    };



    double norm()
    {
        // Calculates the Frobenius norm for the sparse matrix
        int n, N;
        N = sm.size();
        double frob = 0.0;

        for (n = 0; n < N; n++)
            frob += (sm[n].val * sm[n].val);

        frob = sqrt(frob);
        return frob;
    };



    SparseMatrix<T> operator* (double s)
    {
        int n;
        int N = size();
        SparseMatrix<T> B = SparseMatrix<T>(N, nrows, ncols);

        for (n = 0; n < N; n++)
            B.push_back(sm[n].row, sm[n].col, sm[n].val * s);

        B.clean_up();
        return B;
    };


    std::vector<double> operator* (std::vector<double> vec)
    {
        int n, i, j;
        int N = size();
        int P = vec.size();

        if (P != nrows && P != ncols)
            throw VectorSizeError();

        std::vector<double> new_vec = std::vector<double>(P, 0.0);

        for (n = 0; n < N; n++)
            new_vec[sm[n].row] += (sm[n].val * vec[sm[n].col]);

        return new_vec;
    };


    SparseMatrix<T> operator* (SparseMatrix<T> B)
    {
        int i, j;
        int N = size();

        if (B.nrows != nrows && B.ncols != ncols)
            throw VectorSizeError();

        // Put B in col-order
        B.sort_by_col();

        typename std::vector< sparse_matrix_entry<T> >::iterator a_iter, b_iter, a_stor, b_stor, a_beg, b_beg, a_end, b_end;
        a_beg = sm.begin();
        a_end = sm.end();
        a_stor = a_beg;

        b_beg = B.sm.begin();
        b_end = B.sm.end();
        b_stor = b_beg;

        SparseMatrix<T> C = SparseMatrix<T>(N, nrows, ncols);

        double dot_prod;
        for (i = 0; i < C.nrows; i++)
        {
            for (j = 0; j < C.ncols; j++)
            {
                a_iter = find_row(a_stor, a_end, i);
                b_iter = find_col(b_stor, b_end, j);

                a_stor = a_iter;
                b_stor = b_iter;

                if ((a_iter == a_end) || (b_iter == b_end))
                    continue;

                dot_prod = 0.0;

                while(((*a_iter).row == i) && ((*b_iter).col == j) && (a_iter != a_end) && (b_iter != b_end))
                {
                    if ((*a_iter).col == (*b_iter).row)
                    {
                        dot_prod += ((*a_iter).val * (*b_iter).val);
                        a_iter++;
                        b_iter++;
                    }
                    else if ((*a_iter).col < (*b_iter).row)
                        a_iter++;
                    else if ((*b_iter).row < (*a_iter).col)
                        b_iter++;
                    else
                        break;
                };

                // Add entry, if not sparse
                if (dot_prod != 0.0)
                    C.push_back(i, j, dot_prod);

                if ((a_iter == a_beg) || (a_iter == a_end))
                    a_stor = a_beg;

                if ((b_iter == b_beg) || (b_iter == b_end))
                    b_stor = b_beg;
            };
        };

        // Put B back in the right order
        B.sort_by_row();

        C.clean_up();
        return C;
    };



    SparseMatrix<T> operator+ (SparseMatrix<T> B)
    {
        int i, j;
        int N = size();

        if (B.nrows != nrows && B.ncols != ncols)
            throw VectorSizeError();

        typename std::vector< sparse_matrix_entry<T> >::iterator a_iter, b_iter, a_end, b_end;
        a_end = sm.end();
        b_end = B.sm.end();

        SparseMatrix<T> C = SparseMatrix<T>(N + B.size(), nrows, ncols);

        double tmp_sum;
        for (i = 0; i < C.nrows; i++)
        {
            a_iter = find_row(sm.begin(), a_end, i);
            b_iter = find_row(B.sm.begin(), b_end, i);

            // Cover the case where there are no a- or b-entries for row == i
            if ((a_iter == a_end) && (b_iter == b_end))
                continue;

            // Cover the case where there are no b-entries for row == i
            if ((a_iter != a_end) && (b_iter == b_end))
            {
                while((*a_iter).row == i)
                {
                    C.push_back(i, (*a_iter).col, (*a_iter).val);
                    a_iter++;
                };
                continue;
            };

            // Cover the case where there are no a-entries for row == i
            if ((a_iter == a_end) && (b_iter != b_end))
            {
                while((*b_iter).row == i)
                {
                    C.push_back(i, (*b_iter).col, (*b_iter).val);
                    b_iter++;
                };
                continue;
            };

            // cover the case when there are both a and b entries for row == i
            while(((*a_iter).row == i) || ((*b_iter).row == i))
            {
                if (((*a_iter).row == i) && ((*b_iter).row == i))
                {
                    if ((*a_iter).col == (*b_iter).col)
                    {
                        tmp_sum = ((*a_iter).val + (*b_iter).val);
                        if (tmp_sum != 0.0)
                            C.push_back(i, (*a_iter).col, tmp_sum);

                        a_iter++;
                        b_iter++;
                    }
                    else if ((*a_iter).col < (*b_iter).col)
                    {
                        C.push_back(i, (*a_iter).col, (*a_iter).val);
                        a_iter++;
                    }
                    else
                    {
                        C.push_back(i, (*b_iter).col, (*b_iter).val);
                        b_iter++;
                    };
                }
                else if ((*a_iter).row == i)
                {
                    C.push_back(i, (*a_iter).col, (*a_iter).val);
                    a_iter++;
                }
                else if ((*b_iter).row == i)
                {
                    C.push_back(i, (*b_iter).col, (*b_iter).val);
                    b_iter++;
                }
                else
                    break;
            };
        };

        C.clean_up();
        return C;
    };

};



// End bright namespace
};

#endif
