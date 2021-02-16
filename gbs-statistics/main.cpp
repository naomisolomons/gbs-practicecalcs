#include <iostream>
#include <vector>
#include <algorithm>

/* int doublefactorial(int n) // This is from the internet (Rahul Agrawal)
{
    if (n == 0 || n==1)
      return 1;
    return n*doublefactorial(n-2);
}


std::vector<std::vector<int>> recursiveLoop(std::vector<int> indices, int N, std::vector<std::vector<int>> tempmu,
                                            std::vector<std::vector<int>> allmus, int endsize) //Not sure about use of const here
{
    if (tempmu.size() == N)
    {
        allmus.push_back(tempmu);
        tempmu.pop_back();
        tempmu.pop_back();
        tempmu.pop_back();
        for (int i{0}; i<(N/2 - 1); ++i)
        {
            if (allmus.size() == endsize)
                break;
            else if (allmus.size() % doublefactorial(N - 2*i - 1) == 0)
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
        std::vector<int> newlist{};
        for(int i=0; i < indices.size(); i++)
        {
            if (std::find(tempmu.begin(), tempmu.end(), indices[i]) = tempmu.end())
                newlist.push_back(indices[i]);
            else
                ;
        }
        for(int i=1, i<newlist.size(), i++)
        {
            tempmu.push_back(newlist[i]);
            allmus = recursiveLoop(indices, N, tempmu, allmus, endsize)
        }
    }

    return allmus;
}

*/
std::vector<std::vector<int>> perfectMatchingPairs(int N)
{
    std::vector<std::vector<int>> allmus;

    if (N==0)
    {
        allmus = {};
    }
    else if (N==2)
    {
        allmus = {{0,1}};
    }
    else
    {
     /*   int endsize{doublefactorial(N - 1)};
        std::vector<int> listofmodes(N);
        std::iota(listofmodes.begin(), listofmodes.end(), 0);
        allmus = rec_loop(listofmodes, N, std::vector<std::vector<int>> tempmu{},
                            std::vector<std::vector<int>> allmus{}, endsize);
    */

        allmus = {{0,1},{0,2}};  // This is a tester line - needs to be deleted
    }

    return allmus;
}


int main()
{

    // The program is going to produce various statistics associated with Gaussian Boson Sampling.
    // The first thing I am going to do is write the function that produces the perfect matching pairs.

    int numberofmodes{6};
    std::vector<std::vector<int>> PMP{perfectMatchingPairs(numberofmodes)};
    for ( const auto &row : PMP) //This is a method for printing 2D vectors from stackexchange
    {
        for ( const auto &s : row ) std::cout << s << ' ';
        std::cout << std::endl;
    }

    return 0;
}
