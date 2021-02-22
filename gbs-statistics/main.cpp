#include <iostream>
#include <vector>
#include <algorithm>
#include <numeric>
#include <tuple> //This is my way of outputting two variables from the function - is it a good?

std::vector<std::vector<int>> x_matrix{{{0,1},{1,0}}};
//Check to see if next bit is in Oli's code?

//std::tuple<std::vector<std::vector<int>>, std::vector<std::vector<int>>> kappa_new

long long doublefactorial(int n) // This is from the internet (Rahul Agrawal)
// Beware of size issues past 33 - write this into below functions?
{
    if (n == 0 || n==1)
      return 1;
    return n*doublefactorial(n-2);
}

int recursiveLoop(std::vector<std::vector<int>> kappaMatrix, const std::vector<int>& indices, int N, std::vector<int>& tempmu,
                                            int& allmus, int endsize, double& hafnian)
{
    if (tempmu.size() == N)
    {
        // std::cout << tempmu[1];
        // allmus.push_back(tempmu);
        allmus += 1;
        double pairing{1};
        for (int k{0}; k<N/2; ++k)
            pairing *= kappaMatrix[tempmu[2*(k+1) - 2]][tempmu[2*(k+1) - 1]];
        hafnian += pairing;
        tempmu.pop_back();
        tempmu.pop_back();
        tempmu.pop_back();
        for (int i{0}; i<(N/2 - 1); ++i)
        {
            if (allmus == endsize)
                break;
            else if (allmus % doublefactorial(N - 2*i - 1) == 0)
                {
                for (int j{0}; j<(N - 2*i - 2); ++j)
                    tempmu.pop_back();
                break;
                }
            else
                ; // In python we have pass here - a 'do nothing' line
        }
    }
    else
    {
        //std::cout << "We get here";
        std::vector<int> newlist{};
        for(int i=0; i < indices.size(); i++)
        {
            if (std::find(tempmu.begin(), tempmu.end(), indices[i]) == tempmu.end())
                newlist.push_back(indices[i]);
            else
                ;
        }
        tempmu.push_back(newlist[0]);
        for(int i=1; i<newlist.size(); i++)
        {
            tempmu.push_back(newlist[i]);
            //for ( const auto &s : tempmu ) std::cout << s << ' ';
            //std::cout << std::endl;
            recursiveLoop(kappaMatrix, indices, N, tempmu, allmus, endsize, hafnian);
        }
    }

    return 0;
}

double Hafnian(std::vector<std::vector<int>> kappaMatrix) //I haven't been able to check this works yet
{
    auto N{kappaMatrix.size()};
    int allmus{0};
    double hafnian{0};

    if (N==0)
    {
        hafnian = 0;
    }
    else if (N==2)
    {
        //allmus = {{0,1}};
        //Got to fill in this case
    }
    else
    {
        long long endsize{doublefactorial(N - 1)};
        std::vector<int> listofmodes(N);
        std::iota(listofmodes.begin(), listofmodes.end(), 0);
        std::vector<int> tempmu{};
        recursiveLoop(kappaMatrix, listofmodes, N, tempmu, allmus, endsize, hafnian);
    }

    return hafnian;
}



int clickProbability()
{

    return 0;

}

int main()
{
    // The program is going to produce various statistics associated with Gaussian Boson Sampling.
    // The first thing I am going to do is write the function that produces the perfect matching pairs.
    std::vector<std::vector<int>> kappaMatrix{{{0,1,1,0},{1,0,0,1},{1,0,1,0},{0,1,0,1}}};
    double PMP{Hafnian(kappaMatrix)};
    //for ( const auto &row : PMP) //This is a method for printing 2D vectors from stackexchange
    //{
    //    for ( const auto &s : row ) std::cout << s << ' ';
    //   std::cout << std::endl;
    //}
    return 0;
}
