#include "common.h"
#include "matrix.hpp"
#include "prime.h"
#include "rational.h"
using namespace RationalNS;
using namespace PRIME;

std::mutex mutex_output;

class p100
{
// see https://euler.stephan-brumme.com/100/
public:
    long long solve()
    {
      unsigned int tests = 1;
      //std::cin >> tests;

      while (tests--)
      {
        unsigned long long minimum = 1000000000000ULL;

        unsigned long long p = 1;
        unsigned long long q = 2;
       // std::cin >> p >> q >> minimum;

        unsigned long long blue = 0;
        unsigned long long red  = 0;

        // special code for p/q = 1/2
        if (p == 1 && q == 2)
        {
          blue = 15;
          red  =  6;

          // keep going until limit is reached
          while (blue + red < minimum)
          {
            // at first I brute-forced the first solutions and wrote them down
            // then I saw the following relationship for p/q = 1/2:
            //  red(i+1) = 2 * blue(i) + red(i) - 1;
            // blue(i+1) = 2 * red(i+1)
            red   = 2 * blue + red - 1; // seems to be true for most p/q
            blue += 2 * red;            // but this line is not correct for p/q != 1/2
          }
          std::cout << blue << std::endl;
          continue;
        }

        // brute-force smallest solution
        bool found = false;
        for (blue = 2; blue < 100000; blue++)
        {
          // sum = red + blue
          // blue * (blue - 1) / (sum * (sum - 1)) = p / q
          // blue * (blue - 1) * q / p = sum * (sum - 1)
          unsigned long long b2 = blue * (blue - 1);
          b2 *= q;
          // right side must be an integer
          if (b2 % p != 0)
            continue;
          unsigned long long sum2 = b2 / p; // sum2 = sum * (sum - 1)

          // (sum-1)^2 < sum2 < sum^2
          unsigned long long sum  = (unsigned long long)std::sqrt(sum2) + 1;
          // sqrt may have returned a floating-point number
          if (sum * (sum - 1) != sum2)
            continue;

          // now we have the correct solution if minimum is small (< 100000)
          red = sum - blue;
          if (blue + red >= minimum)
          {
            found = true;
            break;
          }
        }

        // failed ? TODO: this means just that my simple search aborted
        if (!found)
        {
          std::cout << "No solution" << std::endl;
          continue;
        }

        // show solution
        std::cout << blue << " " << (red + blue) << std::endl;
        return (long long)blue;
      }

      return 0;
    }
};

class p099
{
//709
public:
    int read(std::vector<long long>& vs)
    {
        // READ
        // 519432,525806
        int cnt=0;

#ifdef _WIN32
    std::ifstream is("./../Data/p099_base_exp.txt", std::ifstream::binary);
#else
    std::ifstream is("./../Data/p099_base_exp.txt", std::ifstream::in);
#endif

        if (is)
        {
            // get length of file:
            is.seekg(0, is.end);
            int length = (int)is.tellg();
            is.seekg(0, is.beg);

            std::string s;
            char* buffer = new char[length];

            is.read(buffer, length);
            if (is)
            {
                int pos = 0;
                while(true)
                {
                    if ((buffer[pos]==',') || (buffer[pos]=='\n') )
                    {
                        if (s.size()>0)
                        {
                            vs.push_back(to_long(s));
                            cnt++;
                            s.clear();
                        }
                    }
                    else //if (buffer[pos]!='"')
                        s=s+buffer[pos];

                    pos += 1;
                    if (pos >= length) break;
                }
                if (s.size()>0)
                {
                    vs.push_back(to_long(s));
                    cnt++;
                }
            }
            else
            {
                std::cout << "error: only " << is.gcount() << " could be read of " << length << std::endl;
                std::cout << "BAD FILE" << std::endl;
                return 0;
            }
            is.close();
        }
        else
        {
            std::cout << "FILE not found" << std::endl;
            return 0;
        }
        return cnt;
    }

    int solve()
    {

      std::map<double, unsigned int> data;
      std::vector<long long> vn;
      int n = read(vn);
      std::cout << "read " << n << " numbers" << std::endl;

      for (unsigned int i = 1; i <= 1000; i++)
      {
        unsigned int base, exponent;
        base = (unsigned int)vn[2*i-2];
        exponent = (unsigned int)vn[2*i-1];
        data[exponent * log(base)] = i;
      }
      std::cout << data.rbegin()->second << std::endl;
      return data.rbegin()->second ;
    }
};

class p098
{
// see https://euler.stephan-brumme.com/98/
public:

    int read(std::vector<std::string>& vs)
    {
        // READ
        int cnt=0;

#ifdef _WIN32
    std::ifstream is("./../Data/p098_words.txt", std::ifstream::binary);
#else
    std::ifstream is("./../Data/p098_words.txt", std::ifstream::in);
#endif

        if (is)
        {
            // get length of file:
            is.seekg(0, is.end);
            int length = (int)is.tellg();
            is.seekg(0, is.beg);

            std::string s;
            char* buffer = new char[length];

            is.read(buffer, length);
            if (is)
            {
                int pos = 0;

                while(true)
                {
                    if ((buffer[pos]==',') || (buffer[pos]=='\n') || (buffer[pos]=='"') )
                    {
                        if (s.size()>0)
                        {
                            vs.push_back(s);
                            cnt++;
                            s.clear();
                        }
                    }
                    else if (buffer[pos]!='"')  s=s+buffer[pos];

                    pos += 1;
                    if (pos >= length) break;
                }
                if (s.size()>0)
                {
                    vs.push_back(s);
                    cnt++;
                }
            }
            else
            {
                std::cout << "error: only " << is.gcount() << " could be read of " << length << std::endl;
                std::cout << "BAD FILE" << std::endl;
                return 0;
            }
            is.close();
        }
        else
        {
            std::cout << "FILE not found" << std::endl;
        }
        return cnt;
    }

    // count digits, two numbers have the same fingerprint if they are permutations of each other
    unsigned long long fingerprint(unsigned long long x)
    {
      unsigned long long result = 0;
      while (x > 0)
      {
        auto digit = x % 10;
        x /= 10;

        result += 1LL << (4 * digit);
      }
      return result;
    }

    // return biggest square if both a and b can be mapped to anagram squares, else return 0
    unsigned long long match(const std::string& a, const std::string& b, const std::vector<unsigned long long>& squares)
    {
      unsigned long long result = 0;
      // try all combinations
      for (auto i : squares)
        for (auto j : squares)
        {
          // don't match a word with itself
          if (i == j)
            continue;

          // convert squares to strings
          auto replaceA = std::to_string(i);
          auto replaceB = std::to_string(j);

          // 1. replace all digits of squareA by the letters of a
          // 2. at the same time, whenever such a digit can be found in squareB replace it by the same letter

          // [digit] => letter
          std::map<char, char> replaceTable;
          bool valid = true;
          for (size_t k = 0; k < replaceA.size(); k++)
          {
            char original = replaceA[k];
            // no replacement rule found ? => abort
            if (replaceTable.count(original) != 0 &&
                replaceTable[original] != a[k])
              valid = false;

            // replacement successful
            replaceTable[original] = a[k];
          }

          // two digits must not map to the same letter, though
          std::set<char> used;
          for (auto x : replaceTable)
          {
            // already used ?
            if (used.count(x.second) != 0)
              valid = false;
            // mark as used
            used.insert(x.second);
          }

          // any constraint violation ?
          if (!valid)
            continue;

          // using that mapping, can "a" be constructed ?
          std::string aa;
          for (auto x : replaceA)
            aa += replaceTable[x];
          if (aa != a)
            continue;

          // using that mapping, can "b" be constructed ?
          std::string bb;
          for (auto x : replaceB)
            bb += replaceTable[x];
          if (bb != b)
            continue;

          // new bigger square ?
          if (result < i)
            result = i;
          if (result < j)
            result = j;
        }

      return result;
    }

    long long solve()
    {
      // find word anagrams: sort letters of each word
      // [sorted letters] => [list of words]
      std::map<std::string, std::vector<std::string>> anagrams;

      // read all words from STDIN and fill "anagram" container
      std::vector<std::string> vs;
      int n = read(vs);
      std::cout << "read " << n << " words" << std::endl;
      for(auto& word : vs)
      {
        auto sorted = word;
        std::sort(sorted.begin(), sorted.end());
        // add to word anagrams
        anagrams[sorted].push_back(word);
      }

      // find longest anagram
      size_t maxDigits = 0;
      for (auto i : anagrams)
        if (i.second.size() > 1) // at least two words share the same letters ?
          if (maxDigits < i.second.front().size())
            maxDigits = i.second.front().size();
      // maxDigits will be 9 for the given input ("INTRODUCE", "REDUCTION")

      unsigned long long maxNumber = 1;
      for (size_t i = 0; i < maxDigits; i++)
        maxNumber *= 10;

      // generate all squares
      // for each square, compute its fingerprint
      std::map<unsigned long long, std::vector<unsigned long long>> permutations;
      std::map<unsigned int,       std::set   <unsigned long long>> fingerprintLength;
      // walk through all square numbers (base^2)
      unsigned long long base = 1;
      while (base*base <= maxNumber)
      {
        auto square = base*base;
        auto id     = fingerprint(square);
        permutations[id].push_back(square);

        auto numDigits = log10(square - 1) + 1;
        fingerprintLength[(unsigned int)numDigits].insert(id);

        base++;
      }

      // only process non-unique words (size > 1)
      unsigned long long result = 0;
      for (auto i : anagrams)
      {
        auto pairs = i.second;
        // no other word with the same letters ?
        if (pairs.size() == 1)
          continue;

        // there is a chance that not all words of a permutation group are squares,
        // there only need to be at least two matching words
        auto length = pairs.front().size();
        // compare each word with each other
        for (size_t i = 0; i < pairs.size(); i++)
          for (size_t j = i + 1; j < pairs.size(); j++)
          {
            // extract all relevant squares
            for (auto id : fingerprintLength[(unsigned int)length])
            {
              // and perform the matching process ...
              auto best = match(pairs[i], pairs[j], permutations[id]);
              // bigger square found ?
              if (result < best)
                result = best;
            }
          }
      }

      std::cout << result << std::endl;
      return (long long)result;
    }
};

long long Euler097()
{
    uinteger_t summod = 0;
    uinteger_t mod = 10000000000;
    uinteger_t m(28433);
    uinteger_t prodmod = 1;

    prodmod *= power_modulo(2, 7830457, mod) * m%mod;
    prodmod = prodmod %mod;

    summod += prodmod + (1%mod);
    summod = summod %mod;
    return (long long)summod;
}


std::map<long long, long long> Tinverse(int power, long long t, long double& tminlog);
long long Euler827(long long N=18)
{
//    power: 1 tminlog:1.15118e-4944 2-4 3-1  modulo: 48
//    power: 2 tminlog:3.8712 2-34 3-1  modulo: 399558677
//    power: 3 tminlog:24.6656 2-10 3-2 5-3 7-1  modulo: 8064000
//    power: 4 tminlog:15.9029 2-7 3-3 5-36 7-1  modulo: 236068631
//    power: 5 tminlog:68.0335 2-205 3-81 7-1  modulo: 166932990
//    power: 6 tminlog:233.029 2-9901 5-50  modulo: 404103978
//    power: 7 tminlog:6943.32 2-330 3-44 5-5 7-15  modulo: 20393780
//    power: 8 tminlog:314.313 2-1674 3-18 5-8 7-9 11-2  modulo: 167433903
//    power: 9 tminlog:1215.29 2-66 3-51 5-5 7-5 11-3 13-3 19-2 23-2  modulo: 106205849
//    power: 10 tminlog:146.602 2-1855 3-131 5-50 7-14 11-3  modulo: 272277415
//    power: 11 tminlog:1544.61 2-91 3-48 5-2046 13-11 17-5 29-5  modulo: 181679380
//    power: 12 tminlog:3467.94 2-411 3-228 5-68 7-20 11-6 13-36  modulo: 89637498
//    power: 13 tminlog:790.453 2-2551 3-1190 5-429 7-35 11-1 19-1 23-1  modulo: 110078609
//    power: 14 tminlog:3842.6 2-3331 3-233 5-50 7-135 11-1 13-14 19-1 23-1 31-1  modulo: 181120505
//    power: 15 tminlog:2955.84 2-2262 3-36 5-120 7-8 11-3 13-105 17-6 29-5 37-3  modulo: 90084057
//    power: 16 tminlog:2137.33 2-484 3-239 5-320 7-8 11-2 13-224 17-176 19-2  modulo: 251556746
//    power: 17 tminlog:2212.52 2-6382 3-243 5-2006 7-168 11-3 13-51 17-5 19-1  modulo: 206982739
//    power: 18 tminlog:8401.19 2-1275268146 3-130691468 7-1  modulo: 368953911
//    SUM modulo 397289979

    std::map<long long, long long> rm;
    uinteger_t prodmod = 1;
    uinteger_t mod = 409120391;
    uinteger_t sum = 0;
    long double tminlog = 0;
    std::stringstream ss;

    ss << std::endl;
    for (size_t j=0;j<(size_t)N;j++)
    {
        ss  << "power: " << j+1 << " tminlog:" << tminlog << " ";
        rm = Tinverse((int) j+1, (long long) pow(10, j+1) , tminlog);

        prodmod = 1;
        for(auto& [a,b] : rm)
        {
            //a exp b
            std::cout << a <<"-"<< b << " ";
            ss << a <<"-"<< b << " ";
            prodmod *= power_modulo(a, b, mod);
        }
        prodmod = prodmod % mod;
        ss << " modulo: "<< prodmod << " ";
        sum += (prodmod % mod);

        ss << std::endl;
        std::cout << std::endl;
        std::cout << std::endl;
    }

    uinteger_t summod = sum % mod;
    std::cout << "SUM modulo " << summod << std::endl;
    ss << "SUM modulo " << summod << std::endl;
    to_file("Euler827\n", ss.str());
    return (long long)summod;
}

long long largest_prime_divisor_4kplus1(long long k)
{
    long long m=1;
    long long t = k;
    for (long long i = 2; i <= t; i++)
    {
        if (t == 1) break;
        if (is_prime(i))
        {
            while(t % i == 0)
            {
                //r.push_back(i);
                if (i>2)
                if ((i-1) % 4 == 0)
                {
                    m *= i;
                }
                t = t/i;;
            }
        }
    }
    return m;
}
std::map<uinteger_t, long long> unique_prime_factors_4k1(uinteger_t k, long long lim)
{
    lim = lim;
    std::map<uinteger_t, long long> vmap;
    uinteger_t t = k;
    for (uinteger_t i = 2; i <= t; i++)
    {
        if (t == 1) break;
        if (is_uprime(i))
        {
            while(t % i == 0)
            {
                if (i>2)
                if ((i-1) % 4 == 0)
                {
                    vmap[i]++;
                }
                t = t/i;;
            }
        }
    }
    return vmap;
}
std::map<uinteger_t, long long> unique_prime_factors_4k3(uinteger_t k, long long lim)
{
    lim = lim;
    std::map<uinteger_t, long long> vmap;
    uinteger_t t = k;
    for (uinteger_t i = 2; i <= t; i++)
    {
        if (t == 1) break;
        if (is_uprime(i))
        {
            while(t % i == 0)
            {
                if (i==3)
                {
                    vmap[i]++;
                }
                else
                {
                    if (i>4)
                    if ((i-3) % 4 == 0)
                    {
                        vmap[i]++;
                    }
                }
                t = t/i;;
            }
        }
    }
    return vmap;
}

// a,b,n
long long NumberPythagoreanTriples(long long n)
{
    //https://www.quora.com/In-which-Pythagorean-triplet-is-41-the-smallest-number
    if (n==1) return 0;
    if (n==2) return 0;
    uinteger_t un = n;
    long long m = largest_prime_divisor_4kplus1(n); // m==1 otherwise
    if (n%2==1)
    {
        return (long long)( (ucountDivisors(un*un) + ucountDivisors(m*m))/2 - 1);
    }
    else
    {
        return (long long)( (ucountDivisors(un*un/4) + ucountDivisors(m*m))/2 - 1);
    }
}

long long product_primes(long long k)
{
    long long p = 3;
    long long r = 3; k--;
    while(k > 0)
    {
        p = next_prime(p, 10000000);
        r = r * p;
        k--;
    }
    return r;
}

//https://mathworld.wolfram.com/PythagoreanTriple.html
uinteger_t L(uinteger_t k)
{
    uinteger_t r = 1;
    std::map<uinteger_t, long long> vmap = unique_uprime_factors(k, 10000000);
    if (vmap[2] > 0)
    {
        r *= (2*vmap[2] - 1);
        for(auto& [p, n] : vmap)
        {
            if (p!=2)
                r *= (2*vmap[p] + 1);
        }
    }
    else
    {
        for(auto& [p, n] : vmap)
        {
            if (p!=2)
                r *= (2*vmap[p] + 1);
        }
    }
    return (r-1)/2;
}
uinteger_t sq_func(uinteger_t k)
{
    //https://mathworld.wolfram.com/SumofSquaresFunction.html
    std::map<uinteger_t, long long> vmap3 = unique_prime_factors_4k3(k, 100000000);
    bool has_half = false;
    for(auto& [p, n] : vmap3)
    {
        if (p!=2)
        if (vmap3[p] % 2 == 1)
        {
            has_half = true;
            break;
        }
    }
    if (has_half) return 0;

    std::map<uinteger_t, long long> vmap1 = unique_prime_factors_4k1(k, 100000000);
    long long b = 1;
    for(auto& [p, n] : vmap1)
    {
        b *= (vmap1[p] + 1);
    }
    return 4*b;
}

uinteger_t H(uinteger_t k)
{
    return (sq_func(k*k)-4)/8;
}
uinteger_t T(uinteger_t k)
{
    uinteger_t r = 1;
    return H(k)+L(k);
}

std::mutex mutex_Tinverse;
bool TCheck(int power,  long long t, long long nn, long long kk, long long ll, long long mm, long long kk11,
                        long long kk13,long long kk17,long long kk19,long long kk23,long long kk29)
{
    power = power;
    long long tt = 0;
    long long tt2 = ( (2*ll+1)*(2*kk13+1)*(2*kk17+1)*(2*kk29+1)-1) / 2; //4k+1 = 5, 13 17, 29
    long long tt_beforeN = (2*kk+1)*(2*ll+1)*(2*mm+1)*(2*kk11+1)*(2*kk13+1)*(2*kk17+1)*(2*kk19+1)*(2*kk23+1)*(2*kk29+1);
    if (nn>0)
    {
        tt  = ( ((2*nn-1)*tt_beforeN) - 1 ) / 2;
    }
    else
    {
        tt  = ( tt_beforeN - 1 ) / 2;
    }

    if (t == tt + tt2) return true;
    return false;
}

//https://www.scribd.com/document/255219409/Albert-H-Beiler-Recreations-in-the-theory-of-numbers-the-queen-of-mathematics-entertains-1966-pdf#
std::map<long long, long long> Tinverse(int power, long long t, long double& tminlog)
{
    std::map<long long, long long> rmap;
    std::map<long long, long long> rbestmap;
    long double tlog;
    tminlog = 99999999999;

    std::cout << "Power "<< power << " ";std::cout << std::endl;

    long long T2_2 =  2*(t+1);
    std::cout << "Divisors of "<< T2_2 << " ";std::cout << std::endl;
    std::vector<long long> v = divisors(T2_2);

    long long factor1; long long factor2;
    long long odd_factor;
    long long even_factor;
    for(size_t i = 0; i< v.size()/2;i++)
    {
        factor1 = v[i];
        factor2 = (T2_2 / v[i]);

        rmap.clear();

        if (factor1%2 == 0) even_factor = factor1;
        else if (factor2%2 == 0) even_factor = factor2;
        else even_factor=0;

        if (factor1%2 == 1) odd_factor = factor1;
        else if (factor2%2 == 1) odd_factor = factor2;
        else odd_factor=0;

        std::vector<long long> vexponenta;
        std::vector<long long> vexponentb;
        if (odd_factor>0)
        {
            std::cout << "prime_factors of "<< odd_factor << " ";std::cout << std::endl;
            std::vector<long long> vp = prime_factors(odd_factor);
            for(size_t j = 0; j< vp.size();j++)
            {
                vexponentb.push_back( (vp[j]-1)/2 );
            }
        }

        if (even_factor>0)
        {
            std::cout << "prime_factors of "<< even_factor-1 << " ";std::cout << std::endl;
            std::vector<long long> vp = prime_factors(even_factor-1);
            long long cnt=0;

            if (even_factor==2) vexponenta.push_back( 1 );
            for(size_t j = 0; j< vp.size();j++)
            {
                {
                    if (cnt==0)
                        vexponenta.push_back( (vp[vp.size() - j - 1] + 1)/2 );
                    else
                        vexponenta.push_back( (vp[vp.size() - j - 1]-1)/2 );
                    cnt++;
                }
            }
        }

        long long prime = 2;
        for(size_t j = 0; j< vexponenta.size();j++)
        {
            if (j == 0)
            {
                if (vexponenta[0] > 1)
                    rmap[prime] = vexponenta[0];
            }
            else
            {
                //4k+3
                while( (prime < 3) || ((prime - 3) % 4 != 0))
                {
                    prime = next_prime(prime, prime*prime);
                }
                rmap[prime] = vexponenta[j];
            }
            prime = next_prime(prime, prime*prime);
        }


        prime = 5;
        for(size_t j = 0; j< vexponentb.size();j++)
        {
            {
                //4k+1
                while( ((prime - 1) % 4 != 0) )
                {
                    prime = next_prime(prime, prime*prime);
                }
                rmap[prime] = vexponentb[vexponentb.size() - j - 1];
                prime = next_prime(prime, prime*prime);
            }
        }

        std::cout << factor1 <<"-"<< factor2 << " : ";

        std::cout << " vb:";
        for(size_t j = 0; j< vexponentb.size();j++)
        {
            std::cout << vexponentb[j] << " ";
        }
        std::cout << " va:";
        for(size_t j = 0; j< vexponenta.size();j++)
        {
            std::cout << vexponenta[j] << " ";
        }

        tlog = 0;
        for(auto& [a,b] : rmap)
        {
            tlog += b * log(a);
            std::cout << a <<"-"<< b << " ";
        }

        if (tlog < tminlog)
        {
            tminlog = tlog;
            rbestmap = rmap;
            std::cout  << "tminlog " << tminlog << " ";
        }
        std::cout << std::endl;
    }
    return rbestmap;
}


uinteger_t Tinverse(int power, long long t, long long& n, long long& k, long long& l, long long& m, long long& k11,
                    long long& k13,long long& k17,long long& k19,long long& k23,long long& k29,
                    long long LIM0, long long LIM, long long LIM2
                    )
{
    long double log2 = std::log(2);
    long double log3 = std::log(3);
    long double log5 = std::log(5);
    long double log7 = std::log(7);
    long double log11= std::log(11);
    long double log13 = std::log(13);
    long double log17 = std::log(17);
    long double log19 = std::log(19);
    long double log23 = std::log(23);
    long double log29 = std::log(29);

    uinteger_t cnt =0;
    uinteger_t r;
    uinteger_t tr;
    uinteger_t rmin = -1;
    long long nmin = 0;
    long long kmin = 0;
    long long lmin = 0;
    long long mmin = 0;
    long long k11min = 0;
    long long k13min = 0;
    long long k17min = 0;
    long long k19min = 0;
    long long k23min = 0;
    long long k29min = 0;
    bool ok = false;
    long double tlog;
    long double tminlog = 99999999999;

    long long tt;
    long long tt_beforeN;
    long long tt2;

    uinteger_t maxprime;

    // CACHE
    std::map<long long, std::vector<uinteger_t>> mapcache;
    std::vector<uinteger_t> vprimes;
    std::map<uinteger_t, bool> mapprocessed;

    // tt2 4k+1 primes 5,13,17,29 => give others
    std::map<long long, std::vector<long long>> map_qp;
    std::map<long long, bool> mapt2ok;
    std::vector<long long> vp;

    for(long long ll=LIM0;ll<=LIM;ll++) // 5 = 4k+1 primes
    {
        // Biggest number bucket of 4k+1 primes since minimize the lowest number found (tminlog)
        {
            const std::lock_guard<std::mutex> lock(mutex_Tinverse);
            std::cout<< "5=>"<<  ll << std::endl;
        }

    for(long long kk13=0;kk13<=LIM2;kk13++)
    {
        if (kk13 > ll) break;


    for(long long kk17=0;kk17<=LIM2;kk17++)
    {
        if (kk17 > ll) break;
        if (kk17 > kk13) break;

    for(long long kk29=0;kk29<=LIM2;kk29++)
    {
        if (kk29 > ll) break;
        if (kk29 > kk13) break;
        if (kk29 > kk17) break;
        if (ok)
        {
            tlog = ll*log5 + kk13*log13 + kk17*log17 + kk29*log29;
            if (tlog > tminlog)
                break;
        }

    for(long long z=0;z<1;z++)
    {

            // GIVEN tt2 deduce all others!
            tt2 = ( (2*ll+1)*(2*kk13+1)*(2*kk17+1)*(2*kk29+1) -1) / 2;
            if (t <= tt2)
            {
                continue;
            }

            long long r4k_1 = (2*ll+1)*(2*kk13+1)*(2*kk17+1)*(2*kk29+1);
            long long temp = 2*(t-tt2) + 1;
            long long r0 = temp; // all r4k_1 * r4k_3
            if (r0 % r4k_1 == 0)
            {
                r0 = r0 / r4k_1; //19*5* *3
            }
            // 2*10k+2 = r4k_1 * r4k_3
            // 10 2 3 1 19*5*7*3 -1/2 + 7-1/2
            long long nn=0;
            for(long long zzn=1;zzn<=(r0+1)/2;zzn++) // toooooo many
            {
                if ((2*zzn-1) <= r0)
                if (r0 % (2*zzn-1) == 0) // all odd divisors
                {
                    nn = zzn;

                    long long r2 = r0 / (2*nn-1);

                    long long kk=0;
                    for(long long zz3=0;zz3<=(r2+1)/2;zz3++)
                    {
                        if (zz3 <= nn)
                        if ((2*zz3+1) <= r2)
                        {
                            if (r2 % (2*zz3+1) == 0)
                            {
                                 kk = zz3;
                            }
                            long long r3 = r2 / (2*kk+1);

                            tlog = nn*log2 + ll*log5 + kk13*log13 + kk17*log17 + kk29*log29;
                            if (tlog > tminlog) break;

                            {

                                for(long long mm=0;mm<=(r3+1)/2;mm++)
                                {
                                    if (mm <= kk)
                                    if ((2*mm+1) <= r3)
                                    if (r3 % (2*mm+1) == 0)
                                    {
                                        long long rmm = r3/(2*mm+1) ;
                                        tlog = mm*log7 + kk*log3 + nn*log2 + ll*log5 + kk13*log13 + kk17*log17 + kk29*log29;
                                        if (tlog > tminlog) break;

                                        for(long long kk11=0;kk11<=(rmm+1)/2;kk11++)
                                        {
                                            if (kk11 <= mm)
                                            if ((2*kk11+1) <= rmm)
                                            if (rmm % (2*kk11+1) == 0)
                                            {
                                                long long rkk11 = rmm/(2*kk11+1) ;
                                                tlog = kk11*log11 + mm*log7 + kk*log3 + nn*log2 + ll*log5 + kk13*log13 + kk17*log17 + kk29*log29;
                                                if (tlog > tminlog) break;

                                                for(long long kk19=0;kk19<=(rkk11+1)/2;kk19++)
                                                {
                                                    if (kk19 <= kk11)
                                                    if ((2*kk19+1) <= rkk11)
                                                    if (rkk11 % (2*kk19+1) == 0)
                                                    {
                                                        long long rkk19 = rkk11/(2*kk19+1) ;
                                                        tlog = kk19*log19 + kk11*log11 + mm*log7 + kk*log3 + nn*log2 + ll*log5 + kk13*log13 + kk17*log17 + kk29*log29;
                                                        if (tlog > tminlog) break;

                                                        for(long long kk23=0;kk23<=(rkk19+1)/2;kk23++)
                                                        {
                                                            if (kk23 <= kk19)
                                                            if ((2*kk23+1) <= rkk19)
                                                            if (rkk19 % (2*kk23+1) == 0)
                                                            {
                                                                long long rkk23 = rkk19/(2*kk23+1) ;

                                                                tt_beforeN = (2*kk+1)*(2*ll+1)*(2*mm+1)*(2*kk11+1)*(2*kk13+1)*(2*kk17+1)*(2*kk19+1)*(2*kk23+1)*(2*kk29+1);
                                                                if (nn>0)
                                                                {
                                                                    tt  = ( ((2*nn-1)*tt_beforeN) - 1 ) / 2;
                                                                }
                                                                else
                                                                {
                                                                    tt  = ( tt_beforeN - 1 ) / 2;
                                                                }

                                                                if (t == tt + tt2)
                                                                if (rkk23 == 1)
                                                                {
                                                                    tlog =  nn*log2 + kk*log3 + ll*log5  + mm*log7 + kk11*log11 +
                                                                        kk13*log13 + kk17*log17+ kk19*log19 + kk23*log23 + kk29*log29;

                                                                    {
                                                                        const std::lock_guard<std::mutex> lock(mutex_Tinverse);
                                                                        std::cout
                                                                            << "Test [p:"<<power<<"]" << "[t=" << t << "] "
                                                                            << " Lim=[" << LIM0 << ","<< LIM << ","<< LIM2 << "] "
                                                                            << "[tt:"<< tt << " tt2:"  << tt2 << "] "
                                                                            << "[5=" << ll
                                                                            << " 13=" << kk13
                                                                            << " 17=" << kk17
                                                                            << " 29=" << kk29
                                                                            << "]"
                                                                            << " 2:"
                                                                            << nn << " 3:"
                                                                            << kk << " 5:"
                                                                            << ll << " 7:"
                                                                            << mm << " 11:"
                                                                            << kk11<< " 13:"
                                                                            << kk13 << " 17:"<< kk17 << " "
                                                                            << kk19 << " 23:"<< kk23 << " 29:"<< kk29  << " rmin:"
                                                                            << rmin << " log:" << tlog
                                                                            << " minlog=" << tminlog
                                                                            << " best: "
                                                                            << nmin << " "
                                                                            << kmin << " "
                                                                            << lmin << " "
                                                                            << mmin << " "
                                                                            << k11min << " "
                                                                            << k13min << " "
                                                                            << k17min << " "
                                                                            << k19min << " "
                                                                            << k23min << " "
                                                                            << k29min << " "
                                                                            << " remain `"<< rkk23
                                                                            << std::endl;
                                                                    }

                                                                    if (t == tt + tt2)
                                                                    if (rkk23 == 1)
                                                                    if ((rmin == -1) || (tlog < tminlog))
                                                                    {
                                                                        if (t < 10000000)
                                                                            rmin =  upow(2, nn)    * upow(3, kk)    * upow(5, ll)    * upow(7, mm) *
                                                                                    upow(11, kk11) * upow(13, kk13) * upow(17, kk17) * upow(19, kk19) *
                                                                                    upow(23, kk23) * upow(29, kk29);
                                                                        else
                                                                            rmin = 1234567;
                                                                        if (rmin > 999999999999999) rmin = 1234567;

                                                                        nmin = nn;
                                                                        kmin = kk;
                                                                        lmin = ll;
                                                                        mmin = mm;
                                                                        k11min = kk11;
                                                                        k13min = kk13;
                                                                        k17min = kk17;
                                                                        k19min = kk19;
                                                                        k23min = kk23;
                                                                        k29min = kk29;
                                                                        tminlog = tlog;
                                                                        ok = true;
                                                                        {
                                                                            const std::lock_guard<std::mutex> lock(mutex_Tinverse);
                                                                            std::cout
                                                                                << "NEW [p:"<<power<<"]" << "[t=" << t << "] "
                                                                                << " Lim=[" << LIM0 << ","<< LIM << ","<< LIM2 << "] "
                                                                                << "[tt:"<< tt << " tt2:"  << tt2 << "] "
                                                                                << "[5=" << ll
                                                                                << " 13=" << kk13
                                                                                << " 17=" << kk17
                                                                                << " 29=" << kk29
                                                                                << "]"
                                                                                << " 2:" << nn << " 3:"
                                                                                << kk << " 5:"
                                                                                << ll << " 7:"
                                                                                << mm << " 11:"
                                                                                << kk11<< " 13:"
                                                                                << kk13 << " 17:"<< kk17 << " "
                                                                                << kk19 << " 23:"<< kk23 << " 29:"<< kk29  << " rmin:"
                                                                                << rmin << " log:" << tlog
                                                                                << " minlog=" << tminlog
                                                                                << " best: "
                                                                                << nmin << " "
                                                                                << kmin << " "
                                                                                << lmin << " "
                                                                                << mmin << " "
                                                                                << k11min << " "
                                                                                << k13min << " "
                                                                                << k17min << " "
                                                                                << k19min << " "
                                                                                << k23min << " "
                                                                                << k29min << " "
                                                                                << std::endl;
                                                                        }
                                                                        //break;
                                                                    }
                                                                }
                                                            }
                                                        }
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
    }}
    }}}
    n = nmin;
    k = kmin;
    l = lmin;
    m = mmin;
    k11 = k11min;
    k13 = k13min;
    k17 = k17min;
    k19 = k19min;
    k23 = k23min;
    k29 = k29min;

    uinteger_t mod = 409120391;
    rmin =  power_modulo(2, nmin, mod) *
            power_modulo(3, kmin, mod) *
            power_modulo(5, lmin, mod) *
            power_modulo(7,  mmin, mod) *
            power_modulo(11, k11min, mod) *
            power_modulo(13, k13min, mod) *
            power_modulo(17, k17min, mod) *
            power_modulo(19, k19min, mod) *
            power_modulo(23, k23min, mod) *
            power_modulo(29, k29min, mod)
            ;
    rmin = rmin % mod;

    //    rmin = upow(2, nmin) * upow(3, kmin) * upow(5, lmin) * upow(7, mmin) * upow(11, k11min)
    //            * upow(13, k13min) * upow(17, k17min) * upow(19, k19min) * upow(23, k23min);
    {
        const std::lock_guard<std::mutex> lock(mutex_Tinverse);

        std::cout << "[p:"<<power<<"]" << "[t=" << t << "] "
            << "Exit "
            << " Lim=[" << LIM0 << ","<< LIM << ","<< LIM2 << "] "
            << n << " "
            << k << " "
            << l << " "
            << m << " "
            << k11 << " " << k13 << " "<< k17 << " "<< k19 << " "<< k23<< " "<< k29;
        std::cout << " result modulo " << rmin << " ";
        std::cout << std::endl;
        std::cout << std::endl;
    }
    return rmin;
}

void run_Tinverse(  int power, uinteger_t& result,
                    long long t, long long& n, long long& k, long long& l, long long& m, long long& k11,
                    long long& k13, long long& k17, long long& k19, long long& k23, long long& k29,
                    long long LIM0, long long LIM, long long LIM2)
{
    uinteger_t r = Tinverse(power, t, n, k, l, m, k11, k13, k17, k19, k23, k29, LIM0, LIM, LIM2);
    result  = r;
}

long long Euler827(long FROM, long long N, long long LIM0, long long LIM, long long LIM2, bool do_sieve)
{
    uinteger_t r;
    uinteger_t sum = 0;
    long long nn; long long kk; long long ll; long long mm;long long k11;
    long long k13; long long k17; long long k19; long long k23; long long k29;
//    std::cout << "T 15 "<< T(15) << std::endl;
//    std::cout << "T 48 "<< T(48) << std::endl;
//    std::cout << "T 8064000 "<< T(8064000) << std::endl;

    int nt = std::thread::hardware_concurrency();
    std::vector<std::thread> vt;
    if (nt >= N - FROM + 1)
    {
        std::vector<std::thread> vt;
        std::vector<uinteger_t> vresult(N+1, 0);
        std::vector<long long> vk(N+1, 0);
        std::vector<long long> vn(N+1, 0);
        std::vector<long long> vl(N+1, 0);
        std::vector<long long> vm(N+1, 0);
        std::vector<long long> v11(N+1, 0);
        std::vector<long long> v13(N+1, 0);
        std::vector<long long> v17(N+1, 0);
        std::vector<long long> v19(N+1, 0);
        std::vector<long long> v23(N+1, 0);
        std::vector<long long> v29(N+1, 0);
        nt = (int)( N - FROM + 1);
        for (int i=0;i< nt;i++)
        {
            long long t = (long long) pow(10, FROM+i);
            vt.push_back(   std::thread(run_Tinverse, FROM+i,
                                        std::ref(vresult[i]), t,
                                        std::ref(vn[i]), std::ref(vk[i]), std::ref(vl[i]), std::ref(vm[i]), std::ref(v11[i]),
                                        std::ref(v13[i]), std::ref(v17[i]), std::ref(v19[i]), std::ref(v23[i]), std::ref(v29[i]),
                                        LIM0, LIM, LIM2
                                        )
                        );
        }
        if (do_sieve)
            vt.push_back( std::thread(fill_bitarray_primes, false));
        for (size_t j=0;j<vt.size();j++)  vt[j].join();

        std::cout << "DONE LIM: "<< LIM << std::endl;
        for (int i=0;i< nt;i++)
        {
            sum += (vresult[i] % 409120391);
        }
        uinteger_t rm = sum % 409120391;
        std::cout << "SUM modulo " << sum << std::endl;

        std::stringstream ss;
        for (int i=0;i< nt;i++)
        {
            long double tlog =  vn[i] *log(2)  + vk[i] *log(3)  + vl[i] *log(5)  + vm[i]*log(7) +
                                v11[i]*log(11) + v13[i]*log(13) + v17[i]*log(17) + v19[i]*log(19)
                                + v23[i]*log(23)+ v29[i]*log(29);

            ss   << FROM+i << " Lim=[" << LIM0 << ","<< LIM << ","<< LIM2 << "] "
                 << vn[i]  << " "  << vk[i]  << " "   << vl[i] << " "  << vm[i]  << " " << v11[i] << " "
                 << v13[i] << " "  << v17[i] << " "   << v19[i]<< " "  << v23[i] << " " << v29[i] << " "
                 << (vresult[i] % 409120391) << " "
                 << tlog
                 << std::endl;
        }
        std::cout << ss.str();
        to_file("Euler827\n", ss.str());
        return (long long)sum;
    }
    else
    {
        for(int i=FROM;i<=N;i++)
        {
            long long t = (long long)pow(10, i);
            r = Tinverse(i, t, nn, kk, ll, mm, k11, k13, k17, k19, k23, k29, LIM0, LIM, LIM2);
            std::cout << "Tinv "<< t << " " << nn <<  " " << kk <<  " " << ll <<  " " << mm << std::endl;
            sum += r;
        }
        uinteger_t rm = sum % 409120391;
        std::cout << rm << " ** " << sum << std::endl;
        return (long long)rm;
    }
}

bool solveSudoku(mat::matrix& board, long& guess)
{
    std::vector<bool> number_remaining(9+1, true); number_remaining[0] = false;
    std::vector<bool> v = number_remaining;
    for(int i=0;i<9;i++)
    for(int j=0;j<9;j++)
    {
        if ( board(i,j) != 0) continue;
        v = number_remaining;
        {
            for(int k=0;k<9;k++)
            {
                // Line i
                if (board(i,k)!=0) v[board(i,k)] = false;
                // Column j
                if (board(k,j)!=0) v[board(k,j)] = false;
            }
            // 3x3
            int grid_start_i = 3 * (int)(i/3);
            int grid_start_j = 3 * (int)(j/3);
            for(int ii=0;ii<3;ii++)
            for(int jj=0;jj<3;jj++)
            {
                if (board(grid_start_i+ii,grid_start_j+jj)!=0)
                    v[board(grid_start_i+ii,grid_start_j+jj)] = false;
            }

            // try remaining at ij
            for(int k=1;k<=9;k++)
            {
                if (v[k] == true)
                {
                    board(i,j) = k;
                    guess++;
                    mat::matrix copyboard = board;
                    if (solveSudoku(copyboard, guess) == true)
                    {
                        board = copyboard;
                        return true;
                    }
                }
            }
        }
        return false;
    }
    return true;
}

long long Euler096(long long N)
{
    //24702
    long long sum=0;
    std::string s;
    std::string line[50][9];
    std::ifstream myfile;
    myfile.open ("./../Data/p096_sudoku.txt", std::ios_base::in);
    if (myfile.good())
    {
        for(int k=0;k<N;k++)
        {
            myfile >> s; myfile >> s;
            for(int i=0;i<9;i++)
            {
                myfile >> s;
                line[k][i] = s;
            }
        }
    }
    else
    {
        std::cout << "ERROR opening file"<< std::endl;
    }
    myfile.close();

    mat::matrix board(9,9);
    for(int b=0;b<N;b++)
    {
        for(int l=0;l<9;l++) // line
        {
            s = line[b][l];
            for(int j=0;j<9;j++)
            {
                board(l,j) = s[j]-'0';
                if ((board(l,j) < 0) || (board(l,j) > 9))
                    std::cout << "ERROR CELL "<< l << " " << j << " " << board(l,j)  << std::endl;
//                else
//                    std::cout << "CELL "<< l << " " << j << " " << board(l,j)  << std::endl;
            }
        }
        long guess = 0;
        if (solveSudoku(board, guess))
        {
            sum += 100*board(0,0)+ 10*board(0,1) + board(0,2) ;
        }
    }
    return sum;
}

long long Euler095(long long N)
{
    //14316
    long long r;
    long long k;
    long long cnt = 0;
    long long maxcnt = 0;
    std::vector<long long> vchain;
    std::map<long long, bool> vmapchain;
    std::vector<long long> vmaxchain;

    for(long long i=2; i < N; i++)
    {
        cnt = 1;
        r = i;
        vchain.clear();
        vmapchain.clear();

        vchain.push_back(r);
        vmapchain[r] = true;

        while(true)
        {
            k = r;
            r = sumofFactors(k);
            if (r >= N) break;
            if (r == 0) break;
            if (r == 1) break;

            if (vmapchain.find(r) != vmapchain.end())
            {
                if (r==i)
                {
                    //std::cout << i << "**" << cnt << " "  << vchain.size() << std::endl;
                    if (cnt > maxcnt)
                    {
                        maxcnt = cnt;
                        vmaxchain = vchain;
                    }
                }
                break;
            }
            else
            {
                cnt++;
                vchain.push_back(r);
                vmapchain[r] = true;
            }
        }
    }

    long long rmin = 99999999999;
    if (vmaxchain.size() > 0)
    {
        for(size_t i=0; i < vmaxchain.size(); i++)
        {
            if (vmaxchain[i] < rmin)
            {
                rmin = vmaxchain[i];
            }
        }
    }
    return rmin;
}

bool is_square(long long k)
{
    long long sr = (long long)sqrt(k);
    return (sr * sr == k);
}

long long Euler094(long long N)
 {
    // Heron's Formula
    // 518408346
    long long r;
    long long sum = 0;
    long long s;
    long long c;

    for(long long i=2; 3*i <= 2*N; i++)
    {
        {
            r = 3*i-1;
            if (r%2 == 0)
            {
                s = r/2;
                if (2*s < N)
                {
                    c = i-1;
                    if (is_square(s*(s - c)))
                        sum+=2*s;
                }
            }
        }
        {
            r = 3*i+1;
            if (r%2 == 0)
            {
                s = r/2;
                if (2*s < N)
                {
                    c = i+1;
                    if (is_square(s*(s - c)))
                        sum+=2*s;
                }
            }
        }
        //if (i%10000 == 0) std::cout << i << " " << area2 << " " << sum << std::endl;
        if (3*i-1 >= N + 6) break;
    }
    return (long long)sum;
 }

long long do_oper(long long a,long long b,long long c,long long d, std::vector<int> vop)
{
    RationalNumber r = RationalNumber(a,1);
    int next = 1;
    std::vector<RationalNumber> vn {RationalNumber(a,1),RationalNumber(b,1),RationalNumber(c,1),RationalNumber(d,1)};

    std::string sop[4] = {"+", "-", "*", "/"};
    for(int i=0;i<3;i++)
    {
        if      (vop[i]==0) r = r + vn[next];
        else if (vop[i]==1) r = r - vn[next];
        else if (vop[i]==2) r = r * vn[next];
        else if (vop[i]==3)
        {
            if (vn[next].getM() != 0) r = r/vn[next];
            else return -1;
        }
        next++;
    }
    if (r.getM().getSign())  return -1;
    if (r.getN() != BigIntegerONE) return -1;

    return r.getM().toLongLong();
}

long long do_bin_oper(long long a,long long b,long long c,long long d, std::vector<int> vop)
{
    RationalNumber r = RationalNumber(a,1);
    RationalNumber r2 = RationalNumber(c,1);
    int next = 1;
    int next2 = 3;
    std::vector<RationalNumber> vn {RationalNumber(a,1),RationalNumber(b,1),RationalNumber(c,1),RationalNumber(d,1)};

    std::string sop[4] = {"+", "-", "*", "/"};

    long long i = 0;
    if      (vop[i]==0) r = r + vn[next];
    else if (vop[i]==1) r = r - vn[next];
    else if (vop[i]==2) r = r * vn[next];
    else if (vop[i]==3)
    {
        if (vn[next].getM() != 0) r = r/vn[next];
        else return -1;
    }

    i = 2;  //skip
    if      (vop[i]==0) r2 = r2 + vn[next2];
    else if (vop[i]==1) r2 = r2- vn[next2];
    else if (vop[i]==2) r2 = r2 * vn[next2];
    else if (vop[i]==3)
    {
        if (vn[next2].getM() != 0) r2 = r2/vn[next2];
        else return -1;
    }

    i = 1;
    if      (vop[i]==0) r = r + r2;
    else if (vop[i]==1) r = r - r2;
    else if (vop[i]==2) r = r * r2;
    else if (vop[i]==3)
    {
        if (r2.getM() != 0) r = r/r2;
        else return -1;
    }

    if (r.getM().getSign())  return -1;
    if (r.getN() != BigIntegerONE) return -1;

    return r.getM().toLongLong();
}


long long Euler093(long long N)
{
    //1258
    N=N;
    std::vector<int> voper;
    std::map<std::vector<int>, bool> map_oper;
    std::map<long long, long long> map_result;
    int max_seq_cnt = 0;
    long long r;
    int cnt;
    std::vector<long long> v;
    std::vector<long long> vmax = {0,1,2,3};

    for(int op_1=0;op_1<4;op_1++)
    for(int op_2=0;op_2<4;op_2++)
    for(int op_3=0;op_3<4;op_3++)
    {
        voper.clear();
        voper.push_back(op_1);
        voper.push_back(op_2);
        voper.push_back(op_3);
        map_oper[voper] = true;
    }

    for(long a=0;  a<10;a++)
    for(long b=a+1;b<10;b++)
    for(long c=b+1;c<10;c++)
    for(long d=c+1;d<10;d++)
    {
        map_result.clear();

        //permutation
        v = {a,b,c,d};
        for(int i = 0; i < 4; i++)
        for(int j = 0; j < 4; j++)
        for(int k = 0; k < 4; k++)
        for(int l = 0; l < 4; l++)
        if ((i!=j) && (i!=k) && (i!=l) && (j!=k) && (j!=l) && (k!=l) )
        {
            for(auto& [vop, n] : map_oper)
            {
                r = do_oper(v[i], v[j], v[k], v[l], vop);
                if (r!=-1)
                {
                    map_result[r]++;
                }
            }

            for(auto& [vop, n] : map_oper)
            {
                r = do_bin_oper(v[i], v[j], v[k], v[l], vop);
                if (r!=-1)
                {
                    map_result[r]++;
                }
            }
        }

        cnt = 0;
        for(size_t i=1;i<map_result.size();i++)
        {
            if (map_result.find(i) != map_result.end())
            {
                cnt++;
            }
            else
            {
                break;
            }
        }

        if (cnt > max_seq_cnt)
        {
            max_seq_cnt = cnt;
            std::vector<long long> v = {a,b,c,d};
            vmax = v;
        }
    }

    return 1000*vmax[0]+100*vmax[1]+10*vmax[2]+vmax[3];
}

int Euler828_op(int ioper, int a, int b)
{
    switch (ioper)
    {
        case 0:return a + b;
        case 1:return a - b;
        case 2:return a * b;
        case 3:return (b != 0 && a % b == 0) ? a / b: std::numeric_limits<int>::min();
    }
    return std::numeric_limits<int>::min();
}

bool Euler828_match_target(std::vector<int> vnum, int n, int target)
{
    if (n < 1) return false;
    if (n == 1) return vnum[0] == target;

    // Recursion: Replace all (i, j with i operation j) then search n-1
    std::vector<int> vnewnum(n-1,0);
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            if (i == j) continue;
            int idx=1;
            for (int k = 0; k < n; k++)
            {
                if (k != i && k != j) vnewnum[idx++] = vnum[k];
            }

            for (int iop = 0; iop < 4; iop++)
            {
                vnewnum[0] = Euler828_op(iop, vnum[i], vnum[j]);
                if (vnewnum[0] == std::numeric_limits<int>::min()) continue;

                if (Euler828_match_target(vnewnum, n - 1, target))
                    return true;
            }
        }
    }
    return false;
}

int Euler828_sum(std::vector<int> v, int n)
{
    int sum = 0;
    for (int i = 0; i < n; i++) sum += v[i];
    return sum;
}

long long Euler828_generate_numbers(std::vector<int> vnum, int target, long long& res_sum, bool out=false)
{
    bool ok;
    long long tsum;
    long long minsum = -1;

    // 2
    for (int i = 0; i < 6; i++)
    for (int j = i+1; j < 6; j++)
    {
        if (i!=j)
        {
            std::vector<int> vnumbers(2,0);
            vnumbers[0] = vnum[i];
            vnumbers[1] = vnum[j];
            ok = Euler828_match_target(vnumbers, 2, target);
            if (ok)
            {
                tsum = Euler828_sum(vnumbers, 2);
                if (minsum == -1) minsum = tsum;
                else if (tsum < minsum) minsum=tsum;
            }
        }
    }

    // 3
    for (int i = 0; i < 6; i++)
    for (int j = i+1; j < 6; j++)
    for (int k = j+1; k < 6; k++)
    {
        if ((i!=j)&&(i!=k)&&(j!=k))
        {
            std::vector<int> vnumbers(3,0);
            vnumbers[0] = vnum[i];
            vnumbers[1] = vnum[j];
            vnumbers[2] = vnum[k];
            ok = Euler828_match_target(vnumbers, 3, target);
            if (ok)
            {
                tsum = Euler828_sum(vnumbers, 3);
                if (minsum == -1) minsum=tsum;
                else if (tsum < minsum) minsum=tsum;
            }
        }
    }

    // 4
    for (int i = 0; i < 6; i++)
    for (int j = i+1; j < 6; j++)
    for (int k = j+1; k < 6; k++)
    for (int l = k+1; l < 6; l++)
    {
        if ((i!=j)&&(i!=k)&&(i!=l) &&(j!=k)&&(j!=l) &&(k!=l) )
        {
            std::vector<int> vnumbers(4,0);
            vnumbers[0] = vnum[i];
            vnumbers[1] = vnum[j];
            vnumbers[2] = vnum[k];
            vnumbers[3] = vnum[l];
            ok = Euler828_match_target(vnumbers, 4, target);
            if (ok)
            {
                tsum = Euler828_sum(vnumbers, 4);
                if (minsum == -1)
                    minsum=tsum;
                else if (tsum < minsum)
                    minsum=tsum;
            }
        }
    }

    // 5
    for (int i = 0; i < 6; i++)
    for (int j = i+1; j < 6; j++)
    for (int k = j+1; k < 6; k++)
    for (int l = k+1; l < 6; l++)
    for (int m = l+1; m < 6; m++)
    {
        if ((i!=j)&&(i!=k)&&(i!=l)&&(i!=m)   &&(j!=k)&&(j!=l)&&(j!=m) &&(k!=l)&&(k!=m) &&(l!=m))
        {
            std::vector<int> vnumbers(5,0);
            vnumbers[0] = vnum[i];
            vnumbers[1] = vnum[j];
            vnumbers[2] = vnum[k];
            vnumbers[3] = vnum[l];
            vnumbers[4] = vnum[m];
            ok = Euler828_match_target(vnumbers, 5, target);
            if (ok)
            {
                tsum = Euler828_sum(vnumbers, 5);
                if (minsum == -1)
                    minsum=tsum;
                else if (tsum < minsum)
                    minsum=tsum;
            }
        }
    }

    // 6
    std::vector<int> vnumbers(6,0);
    vnumbers[0] = vnum[0];
    vnumbers[1] = vnum[1];
    vnumbers[2] = vnum[2];
    vnumbers[3] = vnum[3];
    vnumbers[4] = vnum[4];
    vnumbers[5] = vnum[5];
    ok = Euler828_match_target(vnumbers, 6, target);
    if (ok)
    {
        tsum = Euler828_sum(vnumbers, 6);
        if (minsum == -1) minsum=tsum;
        else if (tsum < minsum)
            minsum=tsum;
    }

    if (minsum == -1) minsum=0;
    res_sum = minsum;

    if (out)
    {
        const std::lock_guard<std::mutex> lock(mutex_output);
        std::cout << target << ": " << minsum << std::endl;
    }
    return minsum;
}

long long Euler828()
{
    // READ
    int num[200 * 7];
    int cnt=0;

#ifdef _WIN32
    std::ifstream is("./../Data/p828_number_challenges.txt", std::ifstream::binary);
#else
    std::ifstream is("./../Data/p828_number_challenges.txt", std::ifstream::in);
#endif

    if (is)
    {
        // get length of file:
        is.seekg(0, is.end);
        int length = (int)is.tellg();
        is.seekg(0, is.beg);

        std::string s;
        char* buffer = new char[length];

        is.read(buffer, length);
        if (is)
        {
            int pos = 0;

            while(true)
            {
                if ((buffer[pos]==',') || (buffer[pos]=='\n') || (buffer[pos]==':') )
                {
                    if (s.size()>0)
                    {
                        num[cnt] = (int)to_long(s);
                        cnt++;
                        s.clear();
                    }
                }
                else if (buffer[pos]!=',')  s=s+buffer[pos];

                pos += 1;
                if (pos >= length) break;
            }
            if (s.size()>0)
            {
                num[cnt] = (int)to_long(s);
            }
        }
        else
        {
            std::cout << "error: only " << is.gcount() << " could be read of " << length << std::endl;
            std::cout << "BAD FILE" << std::endl;
            return 0;
        }
        is.close();
    }

    std::vector<long long> vsum(200, 0);
    std::vector<int> v;
    int target;
    uinteger_t totsummod = 0;
    uinteger_t prodmod = 1;
    uinteger_t umod = 1005075251;

    std::vector<std::thread> vt;
    for (int i = 1; i <= 200; i++)
    {
        v.clear();
        target = num[(i-1)*7];
        for (int j = 0; j < 6; j++) v.push_back( num[(i-1)*7 + j + 1] );

        // Multithread
        vt.push_back(std::thread(Euler828_generate_numbers, v, target, std::ref(vsum[i-1]), false) );
    }
    for (size_t j=0;j<vt.size();j++)  vt[j].join();


    for (int i = 1; i <= 200; i++)
    {
        if (vsum[i-1] != 0)
        {
            prodmod = power_modulo(3, i, umod) * (vsum[i-1] % umod);
            prodmod = prodmod % umod;
            totsummod += prodmod;
        }
    }
    totsummod = totsummod % umod;
    std::cout << totsummod << std::endl;
    return (long long)totsummod;
}


long long ssd(long long n)
{
    long long ans = 0;
    while(n)
    {
        long long d = n%10;
        n = n/10;
        ans += d*d;
    }
    return ans;
}

long long Euler092(long long N)
{
    //8581146
    long long n89 = 0;
    for(long long n = 2; n < N; n++)
    {
        long long s = ssd(n);
        while(s != 1 && s != 89)
        {
            s = ssd(s);
        }
        if(s==89){
          n89++;
        }
    }
    return n89;
}

long long Euler091(long long n)
{
    //14234
    long long t = 0;
    long long m = 0;
    for(long long x = 1; x< n+1;x++)
    {
        for(long long y = 1; y< n+1; y++)
        {
            m = gcd(x, y);
            t += std::min(x*m/y, m*(n-y)/x);
        }
    }
    return t*2 + n*n*3;
}

class p090
{
// Adapted From web
public:
    // all sides of a dice
    typedef std::vector<unsigned int> Dice;

    // each cube has 6 different sides, we can choose any 6 out of 10
    const Dice Sides = { 0,1,2,3,4,5,6,7,8,9 };

    const unsigned int Skip = 0;
    const unsigned int Take = 1;
    // if container[x] is Take, then Sides[x] is part of the dice
    const std::vector<unsigned int> Initial = { Skip,Skip,Skip,Skip, Take,Take,Take,Take,Take,Take };
    const std::vector<unsigned int> Unused = { Take };

    long long solve()
    {
        unsigned int dices = 2;
        unsigned int limit = 9; // up to the 9th square (9^2 = 81)
        //std::cin >> limit >> dices;

        // Hackerrank extended the problem to three dices
        const unsigned int AllDices = 3;
        unsigned int maxSquare = 0;

        // generate all square numbers
        std::vector<unsigned short> squares;
        for (unsigned int i = 1; i <= limit; i++)
        {
            auto reduce = i * i;
            maxSquare = reduce;

            std::vector<unsigned int> digits;
            // no matter what, always generate a three-digit square (maybe with some leading zeros)
            for (unsigned int j = 0; j < AllDices; j++)
            {
                auto digit = reduce % 10;
                reduce /= 10;
                // convert all 9s to 6s
                if (digit == 9)
                    digit = 6;

                digits.push_back(digit);
            }

            // digits in ascending order
            std::sort(digits.begin(), digits.end());
            // convert to a fingerprint
            // e.g. 9^2 = 081 (with leading zero)
            // => 018 (sorted)
            // => 18  (sorted, converted to an integer)
            auto sortedSquare = digits[0] * 100 + digits[1] * 10 + digits[2];
            if (std::find(squares.begin(), squares.end(), sortedSquare) == squares.end())
                squares.push_back(sortedSquare);
        }

        // will contain all solutions
        unsigned int valid = 0;
        // all possible label combinations for first dice
        Dice dice1, dice2, dice3;
        auto open = squares;

        auto permutationDice1 = Initial;
        auto permutationDice2 = Initial;
        auto permutationDice3 = Initial;

        do
        {
            dice1.clear();
            for (size_t i = 0; i < permutationDice1.size(); i++)
                if (permutationDice1[i] == Take)
                    dice1.push_back(Sides[i]);

            // if both 6 and 9 are contained, then they behave identical
            //if (permutationDice1[6] == Take && permutationDice1[9] == Take)
              //dice1.pop_back();

            // second dice is "lexicographically" bigger than or equal to first dice
            permutationDice2 = (dices >= 2 ? permutationDice1 : Unused);
            do
            {
                dice2.clear();
                for (size_t i = 0; i < permutationDice2.size(); i++)
                    if (permutationDice2[i] == Take)
                        dice2.push_back(Sides[i]);

                // some digits need to occur at least twice
                if (maxSquare >= 100)
                {
                    // 100 requires two zeros, so we must already have at least one
                    if (std::count(dice1.begin(), dice1.end(), 0) +
                        std::count(dice2.begin(), dice2.end(), 0) < 1)
                        continue;
                }
                if (maxSquare >= 144)
                {
                    // 144 requires two 4s, so we must have at least one by now
                    if (std::count(dice1.begin(), dice1.end(), 4) +
                        std::count(dice2.begin(), dice2.end(), 4) < 1)
                        continue;
                }

                // if less than three dices are requested then add a dummy dice with one side
                permutationDice3 = (dices == 3 ? permutationDice2 : Unused);
                do
                {
                    dice3.clear();
                    for (size_t i = 0; i < permutationDice3.size(); i++)
                        if (permutationDice3[i] == Take)
                            dice3.push_back(Sides[i]);

                    // simple pre-check
                    unsigned int frequency[10] = { 0,0,0,0,0, 0,0,0,0,0 };
                    for (auto x : dice1)
                        frequency[x]++;
                    for (auto x : dice2)
                        frequency[x]++;
                    for (auto x : dice3)
                        frequency[x]++;

                    // for performance optimization only: reject impossible combinations
                    if (frequency[1] < 1)
                        continue;
                    if (maxSquare >= 4 && frequency[4] < 1)
                        continue;
                    if (maxSquare >= 25 && frequency[2] < 1)
                        continue;
                    if (maxSquare >= 25 && frequency[5] < 1)
                        continue;
                    if (maxSquare >= 36 && frequency[3] < 1)
                        continue;
                    if (maxSquare >= 81 && frequency[8] < 1)
                        continue;
                    if (maxSquare >= 100 && frequency[0] < 2)
                        continue;
                    if (maxSquare >= 144 && frequency[4] < 2)
                        continue;

                    std::vector<unsigned int> matches;
                    // build all combinations and remove any squares we encounter along the way
                    for (auto one : dice1)
                    {
                        // 6 is 9 (upside down)
                        if (one == 9)
                            one = 6;
                        for (auto two : dice2)
                        {
                            // 6 is 9 (upside down)
                            if (two == 9)
                                two = 6;

                            for (auto three : dice3)
                            {
                                // 6 is 9 (upside down)
                                if (three == 9)
                                    three = 6;

                                unsigned int current[AllDices] = { one, two, three };

                                // std::sort is much slower for such small containers
                                if (current[0] > current[1])
                                    std::swap(current[0], current[1]);
                                if (current[1] > current[2])
                                    std::swap(current[1], current[2]);
                                if (current[0] > current[1])
                                    std::swap(current[0], current[1]);

                                auto sortedSquare = 100 * current[0] + 10 * current[1] + current[2];
                                // if successful then another square number was matched
                                auto match = std::find(squares.begin(), squares.end(), sortedSquare);
                                // remove it from the list
                                if (match != squares.end())
                                    matches.push_back(sortedSquare);
                            }
                        }
                    }

                    if (matches.size() < squares.size())
                        continue;

                    std::sort(matches.begin(), matches.end());
                    auto last = std::unique(matches.begin(), matches.end());

                    open = squares;
                    for (auto m = matches.begin(); m != last; m++)
                    {
                        auto match = std::find(open.begin(), open.end(), *m);
                        open.erase(match);
                    }

                    // all squares matched ?
                    if (open.empty())
                        valid++;
                } while (std::next_permutation(permutationDice3.begin(), permutationDice3.end()));
            } while (std::next_permutation(permutationDice2.begin(), permutationDice2.end()));
        } while (std::next_permutation(permutationDice1.begin(), permutationDice1.end()));

        //std::cout << valid;
        return valid;
    }
};

long long Euler090(long long n)
{
    // 1217
    n = n;
    p090 c;
    return c.solve();
}

#ifdef _WIN32
#pragma warning(disable:4996)
#endif
class p089
{
// Adapted From web
public:

    int solve()
    {
        FILE* fp;
        char str[100];
        int counter; int sum = 0;
        #ifdef _WIN32
        fp = fopen("p089_roman.txt", "r");
        #else
        fp = fopen("../Euler/p089_roman.txt", "r");
        #endif
        if (fp != NULL)
        {
            while (!feof(fp))
            {
                auto r = fscanf(fp, "%s", str);
                if (r==0) {}
                counter = (int)strlen(str);
                str[counter] = '\0';
                //printf("%s\n",str);
                sum += counter - findRoman(findNumber(str));
                str[0] = '\0';
            }
        }
        else
            printf("there was an error in opening the file...aborting...");

        return sum;
    }

    int findRoman(int n)
    {
        int counter = 0;
        //printf("%d => ",n);
        while (n >= 1000)
        {
            counter++;
            //printf("M");
            n -= 1000;
        }
        if (n >= 900)
        {
            counter += 2;
            //printf("CM");
            n -= 900;
        }
        if (n >= 500)
        {
            counter++;
            //printf("D");
            n -= 500;
        }
        if (n >= 400)
        {
            counter += 2;
            //printf("CD");
            n -= 400;
        }
        if (n < 400)
            while (n >= 100)
            {
                counter++;
                //printf("C");
                n -= 100;
            }
        if (n >= 90)
        {
            counter += 2;
            //printf("XC");
            n -= 90;
        }
        if (n >= 50)
        {
            counter++;
            //printf("L");
            n -= 50;
        }
        if (n >= 40)
        {
            counter += 2;
            //printf("XL");
            n -= 40;
        }
        if (n < 40)
            while (n >= 10)
            {
                counter++;
                //printf("X");
                n -= 10;
            }
        if (n == 9)
        {
            counter += 2;
            //printf("IX");
            n -= 9;
        }
        if (n >= 5)
        {
            counter++;
            //printf("V");
            n -= 5;
            //printf("%d",n);
        }
        if (n == 4)
        {
            counter += 2;
            //printf("IV");
            n -= 4;
        }
        while (n != 0)
        {
            counter++;
            //printf("I");
            n--;
            //printf("%d\n",n);
        }
        //printf("\n");
        return counter;
    }

    int findNumber(char* str)
    {
        int i = 0, n = 0, x;
        for (x = 0; x < 2; x++)
        {
            if (str[i] && str[i + 1] && str[i] == 'C' && str[i + 1] == 'M')
            {
                i += 2;
                n += 900;
            }

            while (str[i] && str[i] == 'M')
            {
                n += 1000;
                i++;
            }
        }

        for (x = 0; x < 2; x++)
        {
            if (str[i] && str[i + 1] && str[i] == 'C' && str[i + 1] == 'D')
            {
                i += 2;
                n += 400;
            }

            while (str[i] && str[i] == 'D')
            {
                n += 500;
                i++;
            }
        }

        for (x = 0; x < 2; x++)
        {
            if (str[i] && str[i + 1] && str[i] == 'X' && str[i + 1] == 'C')
            {
                i += 2;
                n += 90;
            }

            while (str[i] && str[i] == 'C')
            {
                n += 100;
                i++;
            }
        }

        for (x = 0; x < 2; x++)
        {
            if (str[i] && str[i + 1] && str[i] == 'X' && str[i + 1] == 'L')
            {
                i += 2;
                n += 40;
            }

            while (str[i] && str[i] == 'L')
            {
                n += 50;
                i++;
            }
        }

        for (x = 0; x < 2; x++)
        {
            if (str[i] && str[i + 1] && str[i] == 'I' && str[i + 1] == 'X')
            {
                i += 2;
                n += 9;
            }

            while (str[i] && str[i] == 'X')
            {
                n += 10;
                i++;
            }
        }

        for (x = 0; x < 2; x++)
        {
            if (str[i] && str[i + 1] && str[i] == 'I' && str[i + 1] == 'V')
            {
                i += 2;
                n += 4;
            }

            while (str[i] && str[i] == 'V')
            {
                n += 5;
                i++;
            }
        }

        while (str[i] && str[i] == 'I')
        {
            n++;
            i++;
        }

        return n;
    }
};

long long Euler089(long long n)
{
    // 743
    n = n;
    p089 c;
    return c.solve();
}

class p088
{
// Adapted From web
public:
    const int maxK = 12000;
    std::vector<int> n;

    void mps(int number, int sum, int product, int start) {
        int k = number - sum + product;
        if (k < maxK) {
            if (number < n[k]) {
                n[k] = number;
            }
            for (int i = start; i < maxK / number * 2; ++i) {
                mps(number * i, sum + i, product + 1, i);
            }
        }
    }

    std::string solve() {
        for (int i = 0; i < maxK; ++i) {
            n.push_back(maxK * 2);
        }
        mps(1, 1, 1, 2);
        std::set<int> ns(n.begin() + 2, n.end());
        return std::to_string(std::accumulate(ns.begin(), ns.end(), 0));
    }
};

std::string Euler088(long long n)
{
    // 7587457
    n = n;
    p088 c;
    return c.solve();
}

class p087
{
// Adapted From web
public:
    constexpr static long long fiftyMillion = 50000001;
    constexpr static long long sz = 10001;
    constexpr static long long mx = 101;
    bool p[sz+1];
    bool* ans = nullptr;// [fiftyMillion + 2] ;
    long long primeTable[1230]; long long  ind = 0;

    ~p087()
    {
        if (ans != nullptr) delete[]ans;
    }

    void sieve() {
        ans = new bool[fiftyMillion + 2];
        long i, j;

        for (i = 0; i < fiftyMillion + 2; i++) ans[i] = false;
        for (i = 0; i <= sz; i++) p[i] = false;
        for (i = 4; i < sz; i += 2)
            p[i] = true;

        primeTable[ind++] = 2;
        for (i = 3; i <= mx; i += 2) {
            if (!p[i]) {
                primeTable[ind++] = i;
                for (j = i + i; j <= sz; j += i)
                    p[j] = true;
            }
        }

        for (i = mx + 2; i <= sz; i += 2) {
            if (!p[i])
                primeTable[ind++] = i;
        }
        //printf("Sieve done\n");
    }

    long long solve()
    {
        long long i; long long j; long long k; long long cnt = 0;
        long long sum; long long S; long long Q; long long F;

        sieve();

        for (i = 0; i < ind; i++)
        {
            S = primeTable[i] * primeTable[i];
            if (S > fiftyMillion)
                break;

            for (j = 0; j < ind; j++)
            {
                Q = primeTable[j] * primeTable[j] * primeTable[j];
                if (Q + S > fiftyMillion)
                    break;

                for (k = 0; k < ind; k++)
                {
                    F = primeTable[k] * primeTable[k] * primeTable[k] * primeTable[k];
                    if (F + Q + S > fiftyMillion)
                        break;

                    sum = S + Q + F;
                    if (!ans[sum])
                    {
                        cnt++;
                        ans[sum] = true;
                    }
                }
            }

        }
        return cnt;
    }
};


long long Euler087(long long n)
{
    // 1097343
    /*
    Roughly we have to generate prime numbers up to 10000.
    There are 1229 prime numbers bellow 10000. After generating these prime
    */
    n = n;
    p087 c;
    return c.solve();
}

long double getLength(long long a, long long b) {
    return sqrt(a*a + b*b);
}
bool is_integer(long double a)
{
    return (a == std::floor(a));
}

long long Euler086(long long n)
{
    // Adapted From web
    // 1818
    // Based on https://www.mathblog.dk/project-euler-86-shortest-path-cuboid/
    long long M = 2;
    long long counter = 0;

    while (counter < n)
    {
        M++;
        for (long long baseHeightWidth = 3; baseHeightWidth <= 2 * M; baseHeightWidth++)
        {
            const long double pathLength = getLength(M, baseHeightWidth);
            if (is_integer(pathLength))
            {
                if (baseHeightWidth <= M) {
                    counter += (long long) std::floor(baseHeightWidth / 2);
                }
                else {
                    counter += 1 + M - (long long)std::floor((baseHeightWidth + 1) / 2);
                }
            }
        }
    }
    return M;
}

long long numberRectangles(long long m, long long n)
{
    return (m + 1) * m * (n + 1) * n / 4;
}

std::string Euler085(long long N)
{
    // Adapted From web
    // 101524
    long long bestDiff = std::numeric_limits<long long>::max();
    long long bestArea = -1;
    long long sqr = (long long)sqrt(N);
    for (long long w = 1; w <= sqr; w++)
    {
        for (long long h = 1; h <= sqr; h++)
        {
            long long diff = std::abs(numberRectangles(w, h) - N);
            if (diff < bestDiff)
            {
                bestDiff = diff;
                bestArea = w * h;
            }
        }
    }
    return std::to_string(bestArea);
}

class p083
{
  // Adapted From web
  // 425185
public:
    int distance[80][80];
    constexpr static int h = 80; constexpr static int w = 80;

#ifdef _WIN32
#pragma warning(push)
#pragma warning(disable:6385)
#endif
    std::string run()
    {
        for (int i = 0; i < 80; i++)
            for (int y = 0; y < 80; y++)
                distance[i][y] = 999999;

        // Bellman�Ford algorithm
        distance[0][0] = GRID[0][0];
        for (int i = 0; i < w * h; i++)
        {
            for (int y = 0; y < h; y++)
            {
                for (int x = 0; x < w; x++)
                {
                    int temp = 999999;
                    temp = std::min(getDistance(x - 1, y), temp);
                    temp = std::min(getDistance(x + 1, y), temp);
                    temp = std::min(getDistance(x, y - 1), temp);
                    temp = std::min(getDistance(x, y + 1), temp);
                    distance[y][x] = std::min(GRID[y][x] + temp, distance[y][x]);
                }
            }
        }
        return std::to_string(distance[h - 1][w - 1]);
    }
#ifdef _WIN32
#pragma warning(pop)
#endif

    int getDistance(int x, int y) {
        if (y < 0 || y >= 80 || x < 0 || x >= 80) return 999999;
        else return distance[y][x];
    }

    static constexpr int GRID[80][80] = {
        {4445,2697,5115,718,2209,2212,654,4348,3079,6821,7668,3276,8874,4190,3785,2752,9473,7817,9137,496,7338,3434,7152,4355,4552,7917,7827,2460,2350,691,3514,5880,3145,7633,7199,3783,5066,7487,3285,1084,8985,760,872,8609,8051,1134,9536,5750,9716,9371,7619,5617,275,9721,2997,2698,1887,8825,6372,3014,2113,7122,7050,6775,5948,2758,1219,3539,348,7989,2735,9862,1263,8089,6401,9462,3168,2758,3748,5870},
        {1096,20,1318,7586,5167,2642,1443,5741,7621,7030,5526,4244,2348,4641,9827,2448,6918,5883,3737,300,7116,6531,567,5997,3971,6623,820,6148,3287,1874,7981,8424,7672,7575,6797,6717,1078,5008,4051,8795,5820,346,1851,6463,2117,6058,3407,8211,117,4822,1317,4377,4434,5925,8341,4800,1175,4173,690,8978,7470,1295,3799,8724,3509,9849,618,3320,7068,9633,2384,7175,544,6583,1908,9983,481,4187,9353,9377},
        {9607,7385,521,6084,1364,8983,7623,1585,6935,8551,2574,8267,4781,3834,2764,2084,2669,4656,9343,7709,2203,9328,8004,6192,5856,3555,2260,5118,6504,1839,9227,1259,9451,1388,7909,5733,6968,8519,9973,1663,5315,7571,3035,4325,4283,2304,6438,3815,9213,9806,9536,196,5542,6907,2475,1159,5820,9075,9470,2179,9248,1828,4592,9167,3713,4640,47,3637,309,7344,6955,346,378,9044,8635,7466,5036,9515,6385,9230},
        {7206,3114,7760,1094,6150,5182,7358,7387,4497,955,101,1478,7777,6966,7010,8417,6453,4955,3496,107,449,8271,131,2948,6185,784,5937,8001,6104,8282,4165,3642,710,2390,575,715,3089,6964,4217,192,5949,7006,715,3328,1152,66,8044,4319,1735,146,4818,5456,6451,4113,1063,4781,6799,602,1504,6245,6550,1417,1343,2363,3785,5448,4545,9371,5420,5068,4613,4882,4241,5043,7873,8042,8434,3939,9256,2187},
        {3620,8024,577,9997,7377,7682,1314,1158,6282,6310,1896,2509,5436,1732,9480,706,496,101,6232,7375,2207,2306,110,6772,3433,2878,8140,5933,8688,1399,2210,7332,6172,6403,7333,4044,2291,1790,2446,7390,8698,5723,3678,7104,1825,2040,140,3982,4905,4160,2200,5041,2512,1488,2268,1175,7588,8321,8078,7312,977,5257,8465,5068,3453,3096,1651,7906,253,9250,6021,8791,8109,6651,3412,345,4778,5152,4883,7505},
        {1074,5438,9008,2679,5397,5429,2652,3403,770,9188,4248,2493,4361,8327,9587,707,9525,5913,93,1899,328,2876,3604,673,8576,6908,7659,2544,3359,3883,5273,6587,3065,1749,3223,604,9925,6941,2823,8767,7039,3290,3214,1787,7904,3421,7137,9560,8451,2669,9219,6332,1576,5477,6755,8348,4164,4307,2984,4012,6629,1044,2874,6541,4942,903,1404,9125,5160,8836,4345,2581,460,8438,1538,5507,668,3352,2678,6942},
        {4295,1176,5596,1521,3061,9868,7037,7129,8933,6659,5947,5063,3653,9447,9245,2679,767,714,116,8558,163,3927,8779,158,5093,2447,5782,3967,1716,931,7772,8164,1117,9244,5783,7776,3846,8862,6014,2330,6947,1777,3112,6008,3491,1906,5952,314,4602,8994,5919,9214,3995,5026,7688,6809,5003,3128,2509,7477,110,8971,3982,8539,2980,4689,6343,5411,2992,5270,5247,9260,2269,7474,1042,7162,5206,1232,4556,4757},
        {510,3556,5377,1406,5721,4946,2635,7847,4251,8293,8281,6351,4912,287,2870,3380,3948,5322,3840,4738,9563,1906,6298,3234,8959,1562,6297,8835,7861,239,6618,1322,2553,2213,5053,5446,4402,6500,5182,8585,6900,5756,9661,903,5186,7687,5998,7997,8081,8955,4835,6069,2621,1581,732,9564,1082,1853,5442,1342,520,1737,3703,5321,4793,2776,1508,1647,9101,2499,6891,4336,7012,3329,3212,1442,9993,3988,4930,7706},
        {9444,3401,5891,9716,1228,7107,109,3563,2700,6161,5039,4992,2242,8541,7372,2067,1294,3058,1306,320,8881,5756,9326,411,8650,8824,5495,8282,8397,2000,1228,7817,2099,6473,3571,5994,4447,1299,5991,543,7874,2297,1651,101,2093,3463,9189,6872,6118,872,1008,1779,2805,9084,4048,2123,5877,55,3075,1737,9459,4535,6453,3644,108,5982,4437,5213,1340,6967,9943,5815,669,8074,1838,6979,9132,9315,715,5048},
        {3327,4030,7177,6336,9933,5296,2621,4785,2755,4832,2512,2118,2244,4407,2170,499,7532,9742,5051,7687,970,6924,3527,4694,5145,1306,2165,5940,2425,8910,3513,1909,6983,346,6377,4304,9330,7203,6605,3709,3346,970,369,9737,5811,4427,9939,3693,8436,5566,1977,3728,2399,3985,8303,2492,5366,9802,9193,7296,1033,5060,9144,2766,1151,7629,5169,5995,58,7619,7565,4208,1713,6279,3209,4908,9224,7409,1325,8540},
        {6882,1265,1775,3648,4690,959,5837,4520,5394,1378,9485,1360,4018,578,9174,2932,9890,3696,116,1723,1178,9355,7063,1594,1918,8574,7594,7942,1547,6166,7888,354,6932,4651,1010,7759,6905,661,7689,6092,9292,3845,9605,8443,443,8275,5163,7720,7265,6356,7779,1798,1754,5225,6661,1180,8024,5666,88,9153,1840,3508,1193,4445,2648,3538,6243,6375,8107,5902,5423,2520,1122,5015,6113,8859,9370,966,8673,2442},
        {7338,3423,4723,6533,848,8041,7921,8277,4094,5368,7252,8852,9166,2250,2801,6125,8093,5738,4038,9808,7359,9494,601,9116,4946,2702,5573,2921,9862,1462,1269,2410,4171,2709,7508,6241,7522,615,2407,8200,4189,5492,5649,7353,2590,5203,4274,710,7329,9063,956,8371,3722,4253,4785,1194,4828,4717,4548,940,983,2575,4511,2938,1827,2027,2700,1236,841,5760,1680,6260,2373,3851,1841,4968,1172,5179,7175,3509},
        {4420,1327,3560,2376,6260,2988,9537,4064,4829,8872,9598,3228,1792,7118,9962,9336,4368,9189,6857,1829,9863,6287,7303,7769,2707,8257,2391,2009,3975,4993,3068,9835,3427,341,8412,2134,4034,8511,6421,3041,9012,2983,7289,100,1355,7904,9186,6920,5856,2008,6545,8331,3655,5011,839,8041,9255,6524,3862,8788,62,7455,3513,5003,8413,3918,2076,7960,6108,3638,6999,3436,1441,4858,4181,1866,8731,7745,3744,1000},
        {356,8296,8325,1058,1277,4743,3850,2388,6079,6462,2815,5620,8495,5378,75,4324,3441,9870,1113,165,1544,1179,2834,562,6176,2313,6836,8839,2986,9454,5199,6888,1927,5866,8760,320,1792,8296,7898,6121,7241,5886,5814,2815,8336,1576,4314,3109,2572,6011,2086,9061,9403,3947,5487,9731,7281,3159,1819,1334,3181,5844,5114,9898,4634,2531,4412,6430,4262,8482,4546,4555,6804,2607,9421,686,8649,8860,7794,6672},
        {9870,152,1558,4963,8750,4754,6521,6256,8818,5208,5691,9659,8377,9725,5050,5343,2539,6101,1844,9700,7750,8114,5357,3001,8830,4438,199,9545,8496,43,2078,327,9397,106,6090,8181,8646,6414,7499,5450,4850,6273,5014,4131,7639,3913,6571,8534,9703,4391,7618,445,1320,5,1894,6771,7383,9191,4708,9706,6939,7937,8726,9382,5216,3685,2247,9029,8154,1738,9984,2626,9438,4167,6351,5060,29,1218,1239,4785},
        {192,5213,8297,8974,4032,6966,5717,1179,6523,4679,9513,1481,3041,5355,9303,9154,1389,8702,6589,7818,6336,3539,5538,3094,6646,6702,6266,2759,4608,4452,617,9406,8064,6379,444,5602,4950,1810,8391,1536,316,8714,1178,5182,5863,5110,5372,4954,1978,2971,5680,4863,2255,4630,5723,2168,538,1692,1319,7540,440,6430,6266,7712,7385,5702,620,641,3136,7350,1478,3155,2820,9109,6261,1122,4470,14,8493,2095},
        {1046,4301,6082,474,4974,7822,2102,5161,5172,6946,8074,9716,6586,9962,9749,5015,2217,995,5388,4402,7652,6399,6539,1349,8101,3677,1328,9612,7922,2879,231,5887,2655,508,4357,4964,3554,5930,6236,7384,4614,280,3093,9600,2110,7863,2631,6626,6620,68,1311,7198,7561,1768,5139,1431,221,230,2940,968,5283,6517,2146,1646,869,9402,7068,8645,7058,1765,9690,4152,2926,9504,2939,7504,6074,2944,6470,7859},
        {4659,736,4951,9344,1927,6271,8837,8711,3241,6579,7660,5499,5616,3743,5801,4682,9748,8796,779,1833,4549,8138,4026,775,4170,2432,4174,3741,7540,8017,2833,4027,396,811,2871,1150,9809,2719,9199,8504,1224,540,2051,3519,7982,7367,2761,308,3358,6505,2050,4836,5090,7864,805,2566,2409,6876,3361,8622,5572,5895,3280,441,7893,8105,1634,2929,274,3926,7786,6123,8233,9921,2674,5340,1445,203,4585,3837},
        {5759,338,7444,7968,7742,3755,1591,4839,1705,650,7061,2461,9230,9391,9373,2413,1213,431,7801,4994,2380,2703,6161,6878,8331,2538,6093,1275,5065,5062,2839,582,1014,8109,3525,1544,1569,8622,7944,2905,6120,1564,1839,5570,7579,1318,2677,5257,4418,5601,7935,7656,5192,1864,5886,6083,5580,6202,8869,1636,7907,4759,9082,5854,3185,7631,6854,5872,5632,5280,1431,2077,9717,7431,4256,8261,9680,4487,4752,4286},
        {1571,1428,8599,1230,7772,4221,8523,9049,4042,8726,7567,6736,9033,2104,4879,4967,6334,6716,3994,1269,8995,6539,3610,7667,6560,6065,874,848,4597,1711,7161,4811,6734,5723,6356,6026,9183,2586,5636,1092,7779,7923,8747,6887,7505,9909,1792,3233,4526,3176,1508,8043,720,5212,6046,4988,709,5277,8256,3642,1391,5803,1468,2145,3970,6301,7767,2359,8487,9771,8785,7520,856,1605,8972,2402,2386,991,1383,5963},
        {1822,4824,5957,6511,9868,4113,301,9353,6228,2881,2966,6956,9124,9574,9233,1601,7340,973,9396,540,4747,8590,9535,3650,7333,7583,4806,3593,2738,8157,5215,8472,2284,9473,3906,6982,5505,6053,7936,6074,7179,6688,1564,1103,6860,5839,2022,8490,910,7551,7805,881,7024,1855,9448,4790,1274,3672,2810,774,7623,4223,4850,6071,9975,4935,1915,9771,6690,3846,517,463,7624,4511,614,6394,3661,7409,1395,8127},
        {8738,3850,9555,3695,4383,2378,87,6256,6740,7682,9546,4255,6105,2000,1851,4073,8957,9022,6547,5189,2487,303,9602,7833,1628,4163,6678,3144,8589,7096,8913,5823,4890,7679,1212,9294,5884,2972,3012,3359,7794,7428,1579,4350,7246,4301,7779,7790,3294,9547,4367,3549,1958,8237,6758,3497,3250,3456,6318,1663,708,7714,6143,6890,3428,6853,9334,7992,591,6449,9786,1412,8500,722,5468,1371,108,3939,4199,2535},
        {7047,4323,1934,5163,4166,461,3544,2767,6554,203,6098,2265,9078,2075,4644,6641,8412,9183,487,101,7566,5622,1975,5726,2920,5374,7779,5631,3753,3725,2672,3621,4280,1162,5812,345,8173,9785,1525,955,5603,2215,2580,5261,2765,2990,5979,389,3907,2484,1232,5933,5871,3304,1138,1616,5114,9199,5072,7442,7245,6472,4760,6359,9053,7876,2564,9404,3043,9026,2261,3374,4460,7306,2326,966,828,3274,1712,3446},
        {3975,4565,8131,5800,4570,2306,8838,4392,9147,11,3911,7118,9645,4994,2028,6062,5431,2279,8752,2658,7836,994,7316,5336,7185,3289,1898,9689,2331,5737,3403,1124,2679,3241,7748,16,2724,5441,6640,9368,9081,5618,858,4969,17,2103,6035,8043,7475,2181,939,415,1617,8500,8253,2155,7843,7974,7859,1746,6336,3193,2617,8736,4079,6324,6645,8891,9396,5522,6103,1857,8979,3835,2475,1310,7422,610,8345,7615},
        {9248,5397,5686,2988,3446,4359,6634,9141,497,9176,6773,7448,1907,8454,916,1596,2241,1626,1384,2741,3649,5362,8791,7170,2903,2475,5325,6451,924,3328,522,90,4813,9737,9557,691,2388,1383,4021,1609,9206,4707,5200,7107,8104,4333,9860,5013,1224,6959,8527,1877,4545,7772,6268,621,4915,9349,5970,706,9583,3071,4127,780,8231,3017,9114,3836,7503,2383,1977,4870,8035,2379,9704,1037,3992,3642,1016,4303},
        {5093,138,4639,6609,1146,5565,95,7521,9077,2272,974,4388,2465,2650,722,4998,3567,3047,921,2736,7855,173,2065,4238,1048,5,6847,9548,8632,9194,5942,4777,7910,8971,6279,7253,2516,1555,1833,3184,9453,9053,6897,7808,8629,4877,1871,8055,4881,7639,1537,7701,2508,7564,5845,5023,2304,5396,3193,2955,1088,3801,6203,1748,3737,1276,13,4120,7715,8552,3047,2921,106,7508,304,1280,7140,2567,9135,5266},
        {6237,4607,7527,9047,522,7371,4883,2540,5867,6366,5301,1570,421,276,3361,527,6637,4861,2401,7522,5808,9371,5298,2045,5096,5447,7755,5115,7060,8529,4078,1943,1697,1764,5453,7085,960,2405,739,2100,5800,728,9737,5704,5693,1431,8979,6428,673,7540,6,7773,5857,6823,150,5869,8486,684,5816,9626,7451,5579,8260,3397,5322,6920,1879,2127,2884,5478,4977,9016,6165,6292,3062,5671,5968,78,4619,4763},
        {9905,7127,9390,5185,6923,3721,9164,9705,4341,1031,1046,5127,7376,6528,3248,4941,1178,7889,3364,4486,5358,9402,9158,8600,1025,874,1839,1783,309,9030,1843,845,8398,1433,7118,70,8071,2877,3904,8866,6722,4299,10,1929,5897,4188,600,1889,3325,2485,6473,4474,7444,6992,4846,6166,4441,2283,2629,4352,7775,1101,2214,9985,215,8270,9750,2740,8361,7103,5930,8664,9690,8302,9267,344,2077,1372,1880,9550},
        {5825,8517,7769,2405,8204,1060,3603,7025,478,8334,1997,3692,7433,9101,7294,7498,9415,5452,3850,3508,6857,9213,6807,4412,7310,854,5384,686,4978,892,8651,3241,2743,3801,3813,8588,6701,4416,6990,6490,3197,6838,6503,114,8343,5844,8646,8694,65,791,5979,2687,2621,2019,8097,1423,3644,9764,4921,3266,3662,5561,2476,8271,8138,6147,1168,3340,1998,9874,6572,9873,6659,5609,2711,3931,9567,4143,7833,8887},
        {6223,2099,2700,589,4716,8333,1362,5007,2753,2848,4441,8397,7192,8191,4916,9955,6076,3370,6396,6971,3156,248,3911,2488,4930,2458,7183,5455,170,6809,6417,3390,1956,7188,577,7526,2203,968,8164,479,8699,7915,507,6393,4632,1597,7534,3604,618,3280,6061,9793,9238,8347,568,9645,2070,5198,6482,5000,9212,6655,5961,7513,1323,3872,6170,3812,4146,2736,67,3151,5548,2781,9679,7564,5043,8587,1893,4531},
        {5826,3690,6724,2121,9308,6986,8106,6659,2142,1642,7170,2877,5757,6494,8026,6571,8387,9961,6043,9758,9607,6450,8631,8334,7359,5256,8523,2225,7487,1977,9555,8048,5763,2414,4948,4265,2427,8978,8088,8841,9208,9601,5810,9398,8866,9138,4176,5875,7212,3272,6759,5678,7649,4922,5422,1343,8197,3154,3600,687,1028,4579,2084,9467,4492,7262,7296,6538,7657,7134,2077,1505,7332,6890,8964,4879,7603,7400,5973,739},
        {1861,1613,4879,1884,7334,966,2000,7489,2123,4287,1472,3263,4726,9203,1040,4103,6075,6049,330,9253,4062,4268,1635,9960,577,1320,3195,9628,1030,4092,4979,6474,6393,2799,6967,8687,7724,7392,9927,2085,3200,6466,8702,265,7646,8665,7986,7266,4574,6587,612,2724,704,3191,8323,9523,3002,704,5064,3960,8209,2027,2758,8393,4875,4641,9584,6401,7883,7014,768,443,5490,7506,1852,2005,8850,5776,4487,4269},
        {4052,6687,4705,7260,6645,6715,3706,5504,8672,2853,1136,8187,8203,4016,871,1809,1366,4952,9294,5339,6872,2645,6083,7874,3056,5218,7485,8796,7401,3348,2103,426,8572,4163,9171,3176,948,7654,9344,3217,1650,5580,7971,2622,76,2874,880,2034,9929,1546,2659,5811,3754,7096,7436,9694,9960,7415,2164,953,2360,4194,2397,1047,2196,6827,575,784,2675,8821,6802,7972,5996,6699,2134,7577,2887,1412,4349,4380},
        {4629,2234,6240,8132,7592,3181,6389,1214,266,1910,2451,8784,2790,1127,6932,1447,8986,2492,5476,397,889,3027,7641,5083,5776,4022,185,3364,5701,2442,2840,4160,9525,4828,6602,2614,7447,3711,4505,7745,8034,6514,4907,2605,7753,6958,7270,6936,3006,8968,439,2326,4652,3085,3425,9863,5049,5361,8688,297,7580,8777,7916,6687,8683,7141,306,9569,2384,1500,3346,4601,7329,9040,6097,2727,6314,4501,4974,2829},
        {8316,4072,2025,6884,3027,1808,5714,7624,7880,8528,4205,8686,7587,3230,1139,7273,6163,6986,3914,9309,1464,9359,4474,7095,2212,7302,2583,9462,7532,6567,1606,4436,8981,5612,6796,4385,5076,2007,6072,3678,8331,1338,3299,8845,4783,8613,4071,1232,6028,2176,3990,2148,3748,103,9453,538,6745,9110,926,3125,473,5970,8728,7072,9062,1404,1317,5139,9862,6496,6062,3338,464,1600,2532,1088,8232,7739,8274,3873},
        {2341,523,7096,8397,8301,6541,9844,244,4993,2280,7689,4025,4196,5522,7904,6048,2623,9258,2149,9461,6448,8087,7245,1917,8340,7127,8466,5725,6996,3421,5313,512,9164,9837,9794,8369,4185,1488,7210,1524,1016,4620,9435,2478,7765,8035,697,6677,3724,6988,5853,7662,3895,9593,1185,4727,6025,5734,7665,3070,138,8469,6748,6459,561,7935,8646,2378,462,7755,3115,9690,8877,3946,2728,8793,244,6323,8666,4271},
        {6430,2406,8994,56,1267,3826,9443,7079,7579,5232,6691,3435,6718,5698,4144,7028,592,2627,217,734,6194,8156,9118,58,2640,8069,4127,3285,694,3197,3377,4143,4802,3324,8134,6953,7625,3598,3584,4289,7065,3434,2106,7132,5802,7920,9060,7531,3321,1725,1067,3751,444,5503,6785,7937,6365,4803,198,6266,8177,1470,6390,1606,2904,7555,9834,8667,2033,1723,5167,1666,8546,8152,473,4475,6451,7947,3062,3281},
        {2810,3042,7759,1741,2275,2609,7676,8640,4117,1958,7500,8048,1757,3954,9270,1971,4796,2912,660,5511,3553,1012,5757,4525,6084,7198,8352,5775,7726,8591,7710,9589,3122,4392,6856,5016,749,2285,3356,7482,9956,7348,2599,8944,495,3462,3578,551,4543,7207,7169,7796,1247,4278,6916,8176,3742,8385,2310,1345,8692,2667,4568,1770,8319,3585,4920,3890,4928,7343,5385,9772,7947,8786,2056,9266,3454,2807,877,2660},
        {6206,8252,5928,5837,4177,4333,207,7934,5581,9526,8906,1498,8411,2984,5198,5134,2464,8435,8514,8674,3876,599,5327,826,2152,4084,2433,9327,9697,4800,2728,3608,3849,3861,3498,9943,1407,3991,7191,9110,5666,8434,4704,6545,5944,2357,1163,4995,9619,6754,4200,9682,6654,4862,4744,5953,6632,1054,293,9439,8286,2255,696,8709,1533,1844,6441,430,1999,6063,9431,7018,8057,2920,6266,6799,356,3597,4024,6665},
        {3847,6356,8541,7225,2325,2946,5199,469,5450,7508,2197,9915,8284,7983,6341,3276,3321,16,1321,7608,5015,3362,8491,6968,6818,797,156,2575,706,9516,5344,5457,9210,5051,8099,1617,9951,7663,8253,9683,2670,1261,4710,1068,8753,4799,1228,2621,3275,6188,4699,1791,9518,8701,5932,4275,6011,9877,2933,4182,6059,2930,6687,6682,9771,654,9437,3169,8596,1827,5471,8909,2352,123,4394,3208,8756,5513,6917,2056},
        {5458,8173,3138,3290,4570,4892,3317,4251,9699,7973,1163,1935,5477,6648,9614,5655,9592,975,9118,2194,7322,8248,8413,3462,8560,1907,7810,6650,7355,2939,4973,6894,3933,3784,3200,2419,9234,4747,2208,2207,1945,2899,1407,6145,8023,3484,5688,7686,2737,3828,3704,9004,5190,9740,8643,8650,5358,4426,1522,1707,3613,9887,6956,2447,2762,833,1449,9489,2573,1080,4167,3456,6809,2466,227,7125,2759,6250,6472,8089},
        {3266,7025,9756,3914,1265,9116,7723,9788,6805,5493,2092,8688,6592,9173,4431,4028,6007,7131,4446,4815,3648,6701,759,3312,8355,4485,4187,5188,8746,7759,3528,2177,5243,8379,3838,7233,4607,9187,7216,2190,6967,2920,6082,7910,5354,3609,8958,6949,7731,494,8753,8707,1523,4426,3543,7085,647,6771,9847,646,5049,824,8417,5260,2730,5702,2513,9275,4279,2767,8684,1165,9903,4518,55,9682,8963,6005,2102,6523},
        {1998,8731,936,1479,5259,7064,4085,91,7745,7136,3773,3810,730,8255,2705,2653,9790,6807,2342,355,9344,2668,3690,2028,9679,8102,574,4318,6481,9175,5423,8062,2867,9657,7553,3442,3920,7430,3945,7639,3714,3392,2525,4995,4850,2867,7951,9667,486,9506,9888,781,8866,1702,3795,90,356,1483,4200,2131,6969,5931,486,6880,4404,1084,5169,4910,6567,8335,4686,5043,2614,3352,2667,4513,6472,7471,5720,1616},
        {8878,1613,1716,868,1906,2681,564,665,5995,2474,7496,3432,9491,9087,8850,8287,669,823,347,6194,2264,2592,7871,7616,8508,4827,760,2676,4660,4881,7572,3811,9032,939,4384,929,7525,8419,5556,9063,662,8887,7026,8534,3111,1454,2082,7598,5726,6687,9647,7608,73,3014,5063,670,5461,5631,3367,9796,8475,7908,5073,1565,5008,5295,4457,1274,4788,1728,338,600,8415,8535,9351,7750,6887,5845,1741,125},
        {3637,6489,9634,9464,9055,2413,7824,9517,7532,3577,7050,6186,6980,9365,9782,191,870,2497,8498,2218,2757,5420,6468,586,3320,9230,1034,1393,9886,5072,9391,1178,8464,8042,6869,2075,8275,3601,7715,9470,8786,6475,8373,2159,9237,2066,3264,5000,679,355,3069,4073,494,2308,5512,4334,9438,8786,8637,9774,1169,1949,6594,6072,4270,9158,7916,5752,6794,9391,6301,5842,3285,2141,3898,8027,4310,8821,7079,1307},
        {8497,6681,4732,7151,7060,5204,9030,7157,833,5014,8723,3207,9796,9286,4913,119,5118,7650,9335,809,3675,2597,5144,3945,5090,8384,187,4102,1260,2445,2792,4422,8389,9290,50,1765,1521,6921,8586,4368,1565,5727,7855,2003,4834,9897,5911,8630,5070,1330,7692,7557,7980,6028,5805,9090,8265,3019,3802,698,9149,5748,1965,9658,4417,5994,5584,8226,2937,272,5743,1278,5698,8736,2595,6475,5342,6596,1149,6920},
        {8188,8009,9546,6310,8772,2500,9846,6592,6872,3857,1307,8125,7042,1544,6159,2330,643,4604,7899,6848,371,8067,2062,3200,7295,1857,9505,6936,384,2193,2190,301,8535,5503,1462,7380,5114,4824,8833,1763,4974,8711,9262,6698,3999,2645,6937,7747,1128,2933,3556,7943,2885,3122,9105,5447,418,2899,5148,3699,9021,9501,597,4084,175,1621,1,1079,6067,5812,4326,9914,6633,5394,4233,6728,9084,1864,5863,1225},
        {9935,8793,9117,1825,9542,8246,8437,3331,9128,9675,6086,7075,319,1334,7932,3583,7167,4178,1726,7720,695,8277,7887,6359,5912,1719,2780,8529,1359,2013,4498,8072,1129,9998,1147,8804,9405,6255,1619,2165,7491,1,8882,7378,3337,503,5758,4109,3577,985,3200,7615,8058,5032,1080,6410,6873,5496,1466,2412,9885,5904,4406,3605,8770,4361,6205,9193,1537,9959,214,7260,9566,1685,100,4920,7138,9819,5637,976},
        {3466,9854,985,1078,7222,8888,5466,5379,3578,4540,6853,8690,3728,6351,7147,3134,6921,9692,857,3307,4998,2172,5783,3931,9417,2541,6299,13,787,2099,9131,9494,896,8600,1643,8419,7248,2660,2609,8579,91,6663,5506,7675,1947,6165,4286,1972,9645,3805,1663,1456,8853,5705,9889,7489,1107,383,4044,2969,3343,152,7805,4980,9929,5033,1737,9953,7197,9158,4071,1324,473,9676,3984,9680,3606,8160,7384,5432},
        {1005,4512,5186,3953,2164,3372,4097,3247,8697,3022,9896,4101,3871,6791,3219,2742,4630,6967,7829,5991,6134,1197,1414,8923,8787,1394,8852,5019,7768,5147,8004,8825,5062,9625,7988,1110,3992,7984,9966,6516,6251,8270,421,3723,1432,4830,6935,8095,9059,2214,6483,6846,3120,1587,6201,6691,9096,9627,6671,4002,3495,9939,7708,7465,5879,6959,6634,3241,3401,2355,9061,2611,7830,3941,2177,2146,5089,7079,519,6351},
        {7280,8586,4261,2831,7217,3141,9994,9940,5462,2189,4005,6942,9848,5350,8060,6665,7519,4324,7684,657,9453,9296,2944,6843,7499,7847,1728,9681,3906,6353,5529,2822,3355,3897,7724,4257,7489,8672,4356,3983,1948,6892,7415,4153,5893,4190,621,1736,4045,9532,7701,3671,1211,1622,3176,4524,9317,7800,5638,6644,6943,5463,3531,2821,1347,5958,3436,1438,2999,994,850,4131,2616,1549,3465,5946,690,9273,6954,7991},
        {9517,399,3249,2596,7736,2142,1322,968,7350,1614,468,3346,3265,7222,6086,1661,5317,2582,7959,4685,2807,2917,1037,5698,1529,3972,8716,2634,3301,3412,8621,743,8001,4734,888,7744,8092,3671,8941,1487,5658,7099,2781,99,1932,4443,4756,4652,9328,1581,7855,4312,5976,7255,6480,3996,2748,1973,9731,4530,2790,9417,7186,5303,3557,351,7182,9428,1342,9020,7599,1392,8304,2070,9138,7215,2008,9937,1106,7110},
        {7444,769,9688,632,1571,6820,8743,4338,337,3366,3073,1946,8219,104,4210,6986,249,5061,8693,7960,6546,1004,8857,5997,9352,4338,6105,5008,2556,6518,6694,4345,3727,7956,20,3954,8652,4424,9387,2035,8358,5962,5304,5194,8650,8282,1256,1103,2138,6679,1985,3653,2770,2433,4278,615,2863,1715,242,3790,2636,6998,3088,1671,2239,957,5411,4595,6282,2881,9974,2401,875,7574,2987,4587,3147,6766,9885,2965},
        {3287,3016,3619,6818,9073,6120,5423,557,2900,2015,8111,3873,1314,4189,1846,4399,7041,7583,2427,2864,3525,5002,2069,748,1948,6015,2684,438,770,8367,1663,7887,7759,1885,157,7770,4520,4878,3857,1137,3525,3050,6276,5569,7649,904,4533,7843,2199,5648,7628,9075,9441,3600,7231,2388,5640,9096,958,3058,584,5899,8150,1181,9616,1098,8162,6819,8171,1519,1140,7665,8801,2632,1299,9192,707,9955,2710,7314},
        {1772,2963,7578,3541,3095,1488,7026,2634,6015,4633,4370,2762,1650,2174,909,8158,2922,8467,4198,4280,9092,8856,8835,5457,2790,8574,9742,5054,9547,4156,7940,8126,9824,7340,8840,6574,3547,1477,3014,6798,7134,435,9484,9859,3031,4,1502,4133,1738,1807,4825,463,6343,9701,8506,9822,9555,8688,8168,3467,3234,6318,1787,5591,419,6593,7974,8486,9861,6381,6758,194,3061,4315,2863,4665,3789,2201,1492,4416},
        {126,8927,6608,5682,8986,6867,1715,6076,3159,788,3140,4744,830,9253,5812,5021,7616,8534,1546,9590,1101,9012,9821,8132,7857,4086,1069,7491,2988,1579,2442,4321,2149,7642,6108,250,6086,3167,24,9528,7663,2685,1220,9196,1397,5776,1577,1730,5481,977,6115,199,6326,2183,3767,5928,5586,7561,663,8649,9688,949,5913,9160,1870,5764,9887,4477,6703,1413,4995,5494,7131,2192,8969,7138,3997,8697,646,1028},
        {8074,1731,8245,624,4601,8706,155,8891,309,2552,8208,8452,2954,3124,3469,4246,3352,1105,4509,8677,9901,4416,8191,9283,5625,7120,2952,8881,7693,830,4580,8228,9459,8611,4499,1179,4988,1394,550,2336,6089,6872,269,7213,1848,917,6672,4890,656,1478,6536,3165,4743,4990,1176,6211,7207,5284,9730,4738,1549,4986,4942,8645,3698,9429,1439,2175,6549,3058,6513,1574,6988,8333,3406,5245,5431,7140,7085,6407},
        {7845,4694,2530,8249,290,5948,5509,1588,5940,4495,5866,5021,4626,3979,3296,7589,4854,1998,5627,3926,8346,6512,9608,1918,7070,4747,4182,2858,2766,4606,6269,4107,8982,8568,9053,4244,5604,102,2756,727,5887,2566,7922,44,5986,621,1202,374,6988,4130,3627,6744,9443,4568,1398,8679,397,3928,9159,367,2917,6127,5788,3304,8129,911,2669,1463,9749,264,4478,8940,1109,7309,2462,117,4692,7724,225,2312},
        {4164,3637,2000,941,8903,39,3443,7172,1031,3687,4901,8082,4945,4515,7204,9310,9349,9535,9940,218,1788,9245,2237,1541,5670,6538,6047,5553,9807,8101,1925,8714,445,8332,7309,6830,5786,5736,7306,2710,3034,1838,7969,6318,7912,2584,2080,7437,6705,2254,7428,820,782,9861,7596,3842,3631,8063,5240,6666,394,4565,7865,4895,9890,6028,6117,4724,9156,4473,4552,602,470,6191,4927,5387,884,3146,1978,3000},
        {4258,6880,1696,3582,5793,4923,2119,1155,9056,9698,6603,3768,5514,9927,9609,6166,6566,4536,4985,4934,8076,9062,6741,6163,7399,4562,2337,5600,2919,9012,8459,1308,6072,1225,9306,8818,5886,7243,7365,8792,6007,9256,6699,7171,4230,7002,8720,7839,4533,1671,478,7774,1607,2317,5437,4705,7886,4760,6760,7271,3081,2997,3088,7675,6208,3101,6821,6840,122,9633,4900,2067,8546,4549,2091,7188,5605,8599,6758,5229},
        {7854,5243,9155,3556,8812,7047,2202,1541,5993,4600,4760,713,434,7911,7426,7414,8729,322,803,7960,7563,4908,6285,6291,736,3389,9339,4132,8701,7534,5287,3646,592,3065,7582,2592,8755,6068,8597,1982,5782,1894,2900,6236,4039,6569,3037,5837,7698,700,7815,2491,7272,5878,3083,6778,6639,3589,5010,8313,2581,6617,5869,8402,6808,2951,2321,5195,497,2190,6187,1342,1316,4453,7740,4154,2959,1781,1482,8256},
        {7178,2046,4419,744,8312,5356,6855,8839,319,2962,5662,47,6307,8662,68,4813,567,2712,9931,1678,3101,8227,6533,4933,6656,92,5846,4780,6256,6361,4323,9985,1231,2175,7178,3034,9744,6155,9165,7787,5836,9318,7860,9644,8941,6480,9443,8188,5928,161,6979,2352,5628,6991,1198,8067,5867,6620,3778,8426,2994,3122,3124,6335,3918,8897,2655,9670,634,1088,1576,8935,7255,474,8166,7417,9547,2886,5560,3842},
        {6957,3111,26,7530,7143,1295,1744,6057,3009,1854,8098,5405,2234,4874,9447,2620,9303,27,7410,969,40,2966,5648,7596,8637,4238,3143,3679,7187,690,9980,7085,7714,9373,5632,7526,6707,3951,9734,4216,2146,3602,5371,6029,3039,4433,4855,4151,1449,3376,8009,7240,7027,4602,2947,9081,4045,8424,9352,8742,923,2705,4266,3232,2264,6761,363,2651,3383,7770,6730,7856,7340,9679,2158,610,4471,4608,910,6241},
        {4417,6756,1013,8797,658,8809,5032,8703,7541,846,3357,2920,9817,1745,9980,7593,4667,3087,779,3218,6233,5568,4296,2289,2654,7898,5021,9461,5593,8214,9173,4203,2271,7980,2983,5952,9992,8399,3468,1776,3188,9314,1720,6523,2933,621,8685,5483,8986,6163,3444,9539,4320,155,3992,2828,2150,6071,524,2895,5468,8063,1210,3348,9071,4862,483,9017,4097,6186,9815,3610,5048,1644,1003,9865,9332,2145,1944,2213},
        {9284,3803,4920,1927,6706,4344,7383,4786,9890,2010,5228,1224,3158,6967,8580,8990,8883,5213,76,8306,2031,4980,5639,9519,7184,5645,7769,3259,8077,9130,1317,3096,9624,3818,1770,695,2454,947,6029,3474,9938,3527,5696,4760,7724,7738,2848,6442,5767,6845,8323,4131,2859,7595,2500,4815,3660,9130,8580,7016,8231,4391,8369,3444,4069,4021,556,6154,627,2778,1496,4206,6356,8434,8491,3816,8231,3190,5575,1015},
        {3787,7572,1788,6803,5641,6844,1961,4811,8535,9914,9999,1450,8857,738,4662,8569,6679,2225,7839,8618,286,2648,5342,2294,3205,4546,176,8705,3741,6134,8324,8021,7004,5205,7032,6637,9442,5539,5584,4819,5874,5807,8589,6871,9016,983,1758,3786,1519,6241,185,8398,495,3370,9133,3051,4549,9674,7311,9738,3316,9383,2658,2776,9481,7558,619,3943,3324,6491,4933,153,9738,4623,912,3595,7771,7939,1219,4405},
        {2650,3883,4154,5809,315,7756,4430,1788,4451,1631,6461,7230,6017,5751,138,588,5282,2442,9110,9035,6349,2515,1570,6122,4192,4174,3530,1933,4186,4420,4609,5739,4135,2963,6308,1161,8809,8619,2796,3819,6971,8228,4188,1492,909,8048,2328,6772,8467,7671,9068,2226,7579,6422,7056,8042,3296,2272,3006,2196,7320,3238,3490,3102,37,1293,3212,4767,5041,8773,5794,4456,6174,7279,7054,2835,7053,9088,790,6640},
        {3101,1057,7057,3826,6077,1025,2955,1224,1114,6729,5902,4698,6239,7203,9423,1804,4417,6686,1426,6941,8071,1029,4985,9010,6122,6597,1622,1574,3513,1684,7086,5505,3244,411,9638,4150,907,9135,829,981,1707,5359,8781,9751,5,9131,3973,7159,1340,6955,7514,7993,6964,8198,1933,2797,877,3993,4453,8020,9349,8646,2779,8679,2961,3547,3374,3510,1129,3568,2241,2625,9138,5974,8206,7669,7678,1833,8700,4480},
        {4865,9912,8038,8238,782,3095,8199,1127,4501,7280,2112,2487,3626,2790,9432,1475,6312,8277,4827,2218,5806,7132,8752,1468,7471,6386,739,8762,8323,8120,5169,9078,9058,3370,9560,7987,8585,8531,5347,9312,1058,4271,1159,5286,5404,6925,8606,9204,7361,2415,560,586,4002,2644,1927,2824,768,4409,2942,3345,1002,808,4941,6267,7979,5140,8643,7553,9438,7320,4938,2666,4609,2778,8158,6730,3748,3867,1866,7181},
        {171,3771,7134,8927,4778,2913,3326,2004,3089,7853,1378,1729,4777,2706,9578,1360,5693,3036,1851,7248,2403,2273,8536,6501,9216,613,9671,7131,7719,6425,773,717,8803,160,1114,7554,7197,753,4513,4322,8499,4533,2609,4226,8710,6627,644,9666,6260,4870,5744,7385,6542,6203,7703,6130,8944,5589,2262,6803,6381,7414,6888,5123,7320,9392,9061,6780,322,8975,7050,5089,1061,2260,3199,1150,1865,5386,9699,6501},
        {3744,8454,6885,8277,919,1923,4001,6864,7854,5519,2491,6057,8794,9645,1776,5714,9786,9281,7538,6916,3215,395,2501,9618,4835,8846,9708,2813,3303,1794,8309,7176,2206,1602,1838,236,4593,2245,8993,4017,10,8215,6921,5206,4023,5932,6997,7801,262,7640,3107,8275,4938,7822,2425,3223,3886,2105,8700,9526,2088,8662,8034,7004,5710,2124,7164,3574,6630,9980,4242,2901,9471,1491,2117,4562,1130,9086,4117,6698},
        {2810,2280,2331,1170,4554,4071,8387,1215,2274,9848,6738,1604,7281,8805,439,1298,8318,7834,9426,8603,6092,7944,1309,8828,303,3157,4638,4439,9175,1921,4695,7716,1494,1015,1772,5913,1127,1952,1950,8905,4064,9890,385,9357,7945,5035,7082,5369,4093,6546,5187,5637,2041,8946,1758,7111,6566,1027,1049,5148,7224,7248,296,6169,375,1656,7993,2816,3717,4279,4675,1609,3317,42,6201,3100,3144,163,9530,4531},
        {7096,6070,1009,4988,3538,5801,7149,3063,2324,2912,7911,7002,4338,7880,2481,7368,3516,2016,7556,2193,1388,3865,8125,4637,4096,8114,750,3144,1938,7002,9343,4095,1392,4220,3455,6969,9647,1321,9048,1996,1640,6626,1788,314,9578,6630,2813,6626,4981,9908,7024,4355,3201,3521,3864,3303,464,1923,595,9801,3391,8366,8084,9374,1041,8807,9085,1892,9431,8317,9016,9221,8574,9981,9240,5395,2009,6310,2854,9255},
        {8830,3145,2960,9615,8220,6061,3452,2918,6481,9278,2297,3385,6565,7066,7316,5682,107,7646,4466,68,1952,9603,8615,54,7191,791,6833,2560,693,9733,4168,570,9127,9537,1925,8287,5508,4297,8452,8795,6213,7994,2420,4208,524,5915,8602,8330,2651,8547,6156,1812,6271,7991,9407,9804,1553,6866,1128,2119,4691,9711,8315,5879,9935,6900,482,682,4126,1041,428,6247,3720,5882,7526,2582,4327,7725,3503,2631},
        {2738,9323,721,7434,1453,6294,2957,3786,5722,6019,8685,4386,3066,9057,6860,499,5315,3045,5194,7111,3137,9104,941,586,3066,755,4177,8819,7040,5309,3583,3897,4428,7788,4721,7249,6559,7324,825,7311,3760,6064,6070,9672,4882,584,1365,9739,9331,5783,2624,7889,1604,1303,1555,7125,8312,425,8936,3233,7724,1480,403,7440,1784,1754,4721,1569,652,3893,4574,5692,9730,4813,9844,8291,9199,7101,3391,8914},
        {6044,2928,9332,3328,8588,447,3830,1176,3523,2705,8365,6136,5442,9049,5526,8575,8869,9031,7280,706,2794,8814,5767,4241,7696,78,6570,556,5083,1426,4502,3336,9518,2292,1885,3740,3153,9348,9331,8051,2759,5407,9028,7840,9255,831,515,2612,9747,7435,8964,4971,2048,4900,5967,8271,1719,9670,2810,6777,1594,6367,6259,8316,3815,1689,6840,9437,4361,822,9619,3065,83,6344,7486,8657,8228,9635,6932,4864},
        {8478,4777,6334,4678,7476,4963,6735,3096,5860,1405,5127,7269,7793,4738,227,9168,2996,8928,765,733,1276,7677,6258,1528,9558,3329,302,8901,1422,8277,6340,645,9125,8869,5952,141,8141,1816,9635,4025,4184,3093,83,2344,2747,9352,7966,1206,1126,1826,218,7939,2957,2729,810,8752,5247,4174,4038,8884,7899,9567,301,5265,5752,7524,4381,1669,3106,8270,6228,6373,754,2547,4240,2313,5514,3022,1040,9738},
        {2265,8192,1763,1369,8469,8789,4836,52,1212,6690,5257,8918,6723,6319,378,4039,2421,8555,8184,9577,1432,7139,8078,5452,9628,7579,4161,7490,5159,8559,1011,81,478,5840,1964,1334,6875,8670,9900,739,1514,8692,522,9316,6955,1345,8132,2277,3193,9773,3923,4177,2183,1236,6747,6575,4874,6003,6409,8187,745,8776,9440,7543,9825,2582,7381,8147,7236,5185,7564,6125,218,7991,6394,391,7659,7456,5128,5294},
        {2132,8992,8160,5782,4420,3371,3798,5054,552,5631,7546,4716,1332,6486,7892,7441,4370,6231,4579,2121,8615,1145,9391,1524,1385,2400,9437,2454,7896,7467,2928,8400,3299,4025,7458,4703,7206,6358,792,6200,725,4275,4136,7390,5984,4502,7929,5085,8176,4600,119,3568,76,9363,6943,2248,9077,9731,6213,5817,6729,4190,3092,6910,759,2682,8380,1254,9604,3011,9291,5329,9453,9746,2739,6522,3765,5634,1113,5789},
        {5304,5499,564,2801,679,2653,1783,3608,7359,7797,3284,796,3222,437,7185,6135,8571,2778,7488,5746,678,6140,861,7750,803,9859,9918,2425,3734,2698,9005,4864,9818,6743,2475,132,9486,3825,5472,919,292,4411,7213,7699,6435,9019,6769,1388,802,2124,1345,8493,9487,8558,7061,8777,8833,2427,2238,5409,4957,8503,3171,7622,5779,6145,2417,5873,5563,5693,9574,9491,1937,7384,4563,6842,5432,2751,3406,7981},
    };
};

std::string Euler083()
{
    //425185
    p083 c;
    return c.run();
}

class p082
{
// Adapted From web
public:
	int distance[80][80];

	std::string run()
	{
		int h = 80;
		int w = 80;
		for (int x = 0; x < w; x++)
		{
			for (int y = 0; y < h; y++)
				distance[y][x] = GRID[y][x] + std::min(getValue(x - 1, y), getValue(x, y - 1));

			for (int y = h - 1; y >= 0; y--)
				distance[y][x] = std::min(GRID[y][x] + getValue(x, y + 1), distance[y][x]);
		}

		// Minimum of rightmost column
		int mina = 999999;
		for (int y = 0; y < h; y++)
			mina = std::min(distance[y][w - 1], mina);
		return std::to_string(mina);
	}

	int getValue(int x, int y)
	{
		if (x < 0)
			return 0;
		else if (y < 0 || y >= 80 || x >= 80)
			return 999999;
		else
			return distance[y][x];
	}

	static constexpr int GRID[80][80] = {
		{4445,2697,5115,718,2209,2212,654,4348,3079,6821,7668,3276,8874,4190,3785,2752,9473,7817,9137,496,7338,3434,7152,4355,4552,7917,7827,2460,2350,691,3514,5880,3145,7633,7199,3783,5066,7487,3285,1084,8985,760,872,8609,8051,1134,9536,5750,9716,9371,7619,5617,275,9721,2997,2698,1887,8825,6372,3014,2113,7122,7050,6775,5948,2758,1219,3539,348,7989,2735,9862,1263,8089,6401,9462,3168,2758,3748,5870},
		{1096,20,1318,7586,5167,2642,1443,5741,7621,7030,5526,4244,2348,4641,9827,2448,6918,5883,3737,300,7116,6531,567,5997,3971,6623,820,6148,3287,1874,7981,8424,7672,7575,6797,6717,1078,5008,4051,8795,5820,346,1851,6463,2117,6058,3407,8211,117,4822,1317,4377,4434,5925,8341,4800,1175,4173,690,8978,7470,1295,3799,8724,3509,9849,618,3320,7068,9633,2384,7175,544,6583,1908,9983,481,4187,9353,9377},
		{9607,7385,521,6084,1364,8983,7623,1585,6935,8551,2574,8267,4781,3834,2764,2084,2669,4656,9343,7709,2203,9328,8004,6192,5856,3555,2260,5118,6504,1839,9227,1259,9451,1388,7909,5733,6968,8519,9973,1663,5315,7571,3035,4325,4283,2304,6438,3815,9213,9806,9536,196,5542,6907,2475,1159,5820,9075,9470,2179,9248,1828,4592,9167,3713,4640,47,3637,309,7344,6955,346,378,9044,8635,7466,5036,9515,6385,9230},
		{7206,3114,7760,1094,6150,5182,7358,7387,4497,955,101,1478,7777,6966,7010,8417,6453,4955,3496,107,449,8271,131,2948,6185,784,5937,8001,6104,8282,4165,3642,710,2390,575,715,3089,6964,4217,192,5949,7006,715,3328,1152,66,8044,4319,1735,146,4818,5456,6451,4113,1063,4781,6799,602,1504,6245,6550,1417,1343,2363,3785,5448,4545,9371,5420,5068,4613,4882,4241,5043,7873,8042,8434,3939,9256,2187},
		{3620,8024,577,9997,7377,7682,1314,1158,6282,6310,1896,2509,5436,1732,9480,706,496,101,6232,7375,2207,2306,110,6772,3433,2878,8140,5933,8688,1399,2210,7332,6172,6403,7333,4044,2291,1790,2446,7390,8698,5723,3678,7104,1825,2040,140,3982,4905,4160,2200,5041,2512,1488,2268,1175,7588,8321,8078,7312,977,5257,8465,5068,3453,3096,1651,7906,253,9250,6021,8791,8109,6651,3412,345,4778,5152,4883,7505},
		{1074,5438,9008,2679,5397,5429,2652,3403,770,9188,4248,2493,4361,8327,9587,707,9525,5913,93,1899,328,2876,3604,673,8576,6908,7659,2544,3359,3883,5273,6587,3065,1749,3223,604,9925,6941,2823,8767,7039,3290,3214,1787,7904,3421,7137,9560,8451,2669,9219,6332,1576,5477,6755,8348,4164,4307,2984,4012,6629,1044,2874,6541,4942,903,1404,9125,5160,8836,4345,2581,460,8438,1538,5507,668,3352,2678,6942},
		{4295,1176,5596,1521,3061,9868,7037,7129,8933,6659,5947,5063,3653,9447,9245,2679,767,714,116,8558,163,3927,8779,158,5093,2447,5782,3967,1716,931,7772,8164,1117,9244,5783,7776,3846,8862,6014,2330,6947,1777,3112,6008,3491,1906,5952,314,4602,8994,5919,9214,3995,5026,7688,6809,5003,3128,2509,7477,110,8971,3982,8539,2980,4689,6343,5411,2992,5270,5247,9260,2269,7474,1042,7162,5206,1232,4556,4757},
		{510,3556,5377,1406,5721,4946,2635,7847,4251,8293,8281,6351,4912,287,2870,3380,3948,5322,3840,4738,9563,1906,6298,3234,8959,1562,6297,8835,7861,239,6618,1322,2553,2213,5053,5446,4402,6500,5182,8585,6900,5756,9661,903,5186,7687,5998,7997,8081,8955,4835,6069,2621,1581,732,9564,1082,1853,5442,1342,520,1737,3703,5321,4793,2776,1508,1647,9101,2499,6891,4336,7012,3329,3212,1442,9993,3988,4930,7706},
		{9444,3401,5891,9716,1228,7107,109,3563,2700,6161,5039,4992,2242,8541,7372,2067,1294,3058,1306,320,8881,5756,9326,411,8650,8824,5495,8282,8397,2000,1228,7817,2099,6473,3571,5994,4447,1299,5991,543,7874,2297,1651,101,2093,3463,9189,6872,6118,872,1008,1779,2805,9084,4048,2123,5877,55,3075,1737,9459,4535,6453,3644,108,5982,4437,5213,1340,6967,9943,5815,669,8074,1838,6979,9132,9315,715,5048},
		{3327,4030,7177,6336,9933,5296,2621,4785,2755,4832,2512,2118,2244,4407,2170,499,7532,9742,5051,7687,970,6924,3527,4694,5145,1306,2165,5940,2425,8910,3513,1909,6983,346,6377,4304,9330,7203,6605,3709,3346,970,369,9737,5811,4427,9939,3693,8436,5566,1977,3728,2399,3985,8303,2492,5366,9802,9193,7296,1033,5060,9144,2766,1151,7629,5169,5995,58,7619,7565,4208,1713,6279,3209,4908,9224,7409,1325,8540},
		{6882,1265,1775,3648,4690,959,5837,4520,5394,1378,9485,1360,4018,578,9174,2932,9890,3696,116,1723,1178,9355,7063,1594,1918,8574,7594,7942,1547,6166,7888,354,6932,4651,1010,7759,6905,661,7689,6092,9292,3845,9605,8443,443,8275,5163,7720,7265,6356,7779,1798,1754,5225,6661,1180,8024,5666,88,9153,1840,3508,1193,4445,2648,3538,6243,6375,8107,5902,5423,2520,1122,5015,6113,8859,9370,966,8673,2442},
		{7338,3423,4723,6533,848,8041,7921,8277,4094,5368,7252,8852,9166,2250,2801,6125,8093,5738,4038,9808,7359,9494,601,9116,4946,2702,5573,2921,9862,1462,1269,2410,4171,2709,7508,6241,7522,615,2407,8200,4189,5492,5649,7353,2590,5203,4274,710,7329,9063,956,8371,3722,4253,4785,1194,4828,4717,4548,940,983,2575,4511,2938,1827,2027,2700,1236,841,5760,1680,6260,2373,3851,1841,4968,1172,5179,7175,3509},
		{4420,1327,3560,2376,6260,2988,9537,4064,4829,8872,9598,3228,1792,7118,9962,9336,4368,9189,6857,1829,9863,6287,7303,7769,2707,8257,2391,2009,3975,4993,3068,9835,3427,341,8412,2134,4034,8511,6421,3041,9012,2983,7289,100,1355,7904,9186,6920,5856,2008,6545,8331,3655,5011,839,8041,9255,6524,3862,8788,62,7455,3513,5003,8413,3918,2076,7960,6108,3638,6999,3436,1441,4858,4181,1866,8731,7745,3744,1000},
		{356,8296,8325,1058,1277,4743,3850,2388,6079,6462,2815,5620,8495,5378,75,4324,3441,9870,1113,165,1544,1179,2834,562,6176,2313,6836,8839,2986,9454,5199,6888,1927,5866,8760,320,1792,8296,7898,6121,7241,5886,5814,2815,8336,1576,4314,3109,2572,6011,2086,9061,9403,3947,5487,9731,7281,3159,1819,1334,3181,5844,5114,9898,4634,2531,4412,6430,4262,8482,4546,4555,6804,2607,9421,686,8649,8860,7794,6672},
		{9870,152,1558,4963,8750,4754,6521,6256,8818,5208,5691,9659,8377,9725,5050,5343,2539,6101,1844,9700,7750,8114,5357,3001,8830,4438,199,9545,8496,43,2078,327,9397,106,6090,8181,8646,6414,7499,5450,4850,6273,5014,4131,7639,3913,6571,8534,9703,4391,7618,445,1320,5,1894,6771,7383,9191,4708,9706,6939,7937,8726,9382,5216,3685,2247,9029,8154,1738,9984,2626,9438,4167,6351,5060,29,1218,1239,4785},
		{192,5213,8297,8974,4032,6966,5717,1179,6523,4679,9513,1481,3041,5355,9303,9154,1389,8702,6589,7818,6336,3539,5538,3094,6646,6702,6266,2759,4608,4452,617,9406,8064,6379,444,5602,4950,1810,8391,1536,316,8714,1178,5182,5863,5110,5372,4954,1978,2971,5680,4863,2255,4630,5723,2168,538,1692,1319,7540,440,6430,6266,7712,7385,5702,620,641,3136,7350,1478,3155,2820,9109,6261,1122,4470,14,8493,2095},
		{1046,4301,6082,474,4974,7822,2102,5161,5172,6946,8074,9716,6586,9962,9749,5015,2217,995,5388,4402,7652,6399,6539,1349,8101,3677,1328,9612,7922,2879,231,5887,2655,508,4357,4964,3554,5930,6236,7384,4614,280,3093,9600,2110,7863,2631,6626,6620,68,1311,7198,7561,1768,5139,1431,221,230,2940,968,5283,6517,2146,1646,869,9402,7068,8645,7058,1765,9690,4152,2926,9504,2939,7504,6074,2944,6470,7859},
		{4659,736,4951,9344,1927,6271,8837,8711,3241,6579,7660,5499,5616,3743,5801,4682,9748,8796,779,1833,4549,8138,4026,775,4170,2432,4174,3741,7540,8017,2833,4027,396,811,2871,1150,9809,2719,9199,8504,1224,540,2051,3519,7982,7367,2761,308,3358,6505,2050,4836,5090,7864,805,2566,2409,6876,3361,8622,5572,5895,3280,441,7893,8105,1634,2929,274,3926,7786,6123,8233,9921,2674,5340,1445,203,4585,3837},
		{5759,338,7444,7968,7742,3755,1591,4839,1705,650,7061,2461,9230,9391,9373,2413,1213,431,7801,4994,2380,2703,6161,6878,8331,2538,6093,1275,5065,5062,2839,582,1014,8109,3525,1544,1569,8622,7944,2905,6120,1564,1839,5570,7579,1318,2677,5257,4418,5601,7935,7656,5192,1864,5886,6083,5580,6202,8869,1636,7907,4759,9082,5854,3185,7631,6854,5872,5632,5280,1431,2077,9717,7431,4256,8261,9680,4487,4752,4286},
		{1571,1428,8599,1230,7772,4221,8523,9049,4042,8726,7567,6736,9033,2104,4879,4967,6334,6716,3994,1269,8995,6539,3610,7667,6560,6065,874,848,4597,1711,7161,4811,6734,5723,6356,6026,9183,2586,5636,1092,7779,7923,8747,6887,7505,9909,1792,3233,4526,3176,1508,8043,720,5212,6046,4988,709,5277,8256,3642,1391,5803,1468,2145,3970,6301,7767,2359,8487,9771,8785,7520,856,1605,8972,2402,2386,991,1383,5963},
		{1822,4824,5957,6511,9868,4113,301,9353,6228,2881,2966,6956,9124,9574,9233,1601,7340,973,9396,540,4747,8590,9535,3650,7333,7583,4806,3593,2738,8157,5215,8472,2284,9473,3906,6982,5505,6053,7936,6074,7179,6688,1564,1103,6860,5839,2022,8490,910,7551,7805,881,7024,1855,9448,4790,1274,3672,2810,774,7623,4223,4850,6071,9975,4935,1915,9771,6690,3846,517,463,7624,4511,614,6394,3661,7409,1395,8127},
		{8738,3850,9555,3695,4383,2378,87,6256,6740,7682,9546,4255,6105,2000,1851,4073,8957,9022,6547,5189,2487,303,9602,7833,1628,4163,6678,3144,8589,7096,8913,5823,4890,7679,1212,9294,5884,2972,3012,3359,7794,7428,1579,4350,7246,4301,7779,7790,3294,9547,4367,3549,1958,8237,6758,3497,3250,3456,6318,1663,708,7714,6143,6890,3428,6853,9334,7992,591,6449,9786,1412,8500,722,5468,1371,108,3939,4199,2535},
		{7047,4323,1934,5163,4166,461,3544,2767,6554,203,6098,2265,9078,2075,4644,6641,8412,9183,487,101,7566,5622,1975,5726,2920,5374,7779,5631,3753,3725,2672,3621,4280,1162,5812,345,8173,9785,1525,955,5603,2215,2580,5261,2765,2990,5979,389,3907,2484,1232,5933,5871,3304,1138,1616,5114,9199,5072,7442,7245,6472,4760,6359,9053,7876,2564,9404,3043,9026,2261,3374,4460,7306,2326,966,828,3274,1712,3446},
		{3975,4565,8131,5800,4570,2306,8838,4392,9147,11,3911,7118,9645,4994,2028,6062,5431,2279,8752,2658,7836,994,7316,5336,7185,3289,1898,9689,2331,5737,3403,1124,2679,3241,7748,16,2724,5441,6640,9368,9081,5618,858,4969,17,2103,6035,8043,7475,2181,939,415,1617,8500,8253,2155,7843,7974,7859,1746,6336,3193,2617,8736,4079,6324,6645,8891,9396,5522,6103,1857,8979,3835,2475,1310,7422,610,8345,7615},
		{9248,5397,5686,2988,3446,4359,6634,9141,497,9176,6773,7448,1907,8454,916,1596,2241,1626,1384,2741,3649,5362,8791,7170,2903,2475,5325,6451,924,3328,522,90,4813,9737,9557,691,2388,1383,4021,1609,9206,4707,5200,7107,8104,4333,9860,5013,1224,6959,8527,1877,4545,7772,6268,621,4915,9349,5970,706,9583,3071,4127,780,8231,3017,9114,3836,7503,2383,1977,4870,8035,2379,9704,1037,3992,3642,1016,4303},
		{5093,138,4639,6609,1146,5565,95,7521,9077,2272,974,4388,2465,2650,722,4998,3567,3047,921,2736,7855,173,2065,4238,1048,5,6847,9548,8632,9194,5942,4777,7910,8971,6279,7253,2516,1555,1833,3184,9453,9053,6897,7808,8629,4877,1871,8055,4881,7639,1537,7701,2508,7564,5845,5023,2304,5396,3193,2955,1088,3801,6203,1748,3737,1276,13,4120,7715,8552,3047,2921,106,7508,304,1280,7140,2567,9135,5266},
		{6237,4607,7527,9047,522,7371,4883,2540,5867,6366,5301,1570,421,276,3361,527,6637,4861,2401,7522,5808,9371,5298,2045,5096,5447,7755,5115,7060,8529,4078,1943,1697,1764,5453,7085,960,2405,739,2100,5800,728,9737,5704,5693,1431,8979,6428,673,7540,6,7773,5857,6823,150,5869,8486,684,5816,9626,7451,5579,8260,3397,5322,6920,1879,2127,2884,5478,4977,9016,6165,6292,3062,5671,5968,78,4619,4763},
		{9905,7127,9390,5185,6923,3721,9164,9705,4341,1031,1046,5127,7376,6528,3248,4941,1178,7889,3364,4486,5358,9402,9158,8600,1025,874,1839,1783,309,9030,1843,845,8398,1433,7118,70,8071,2877,3904,8866,6722,4299,10,1929,5897,4188,600,1889,3325,2485,6473,4474,7444,6992,4846,6166,4441,2283,2629,4352,7775,1101,2214,9985,215,8270,9750,2740,8361,7103,5930,8664,9690,8302,9267,344,2077,1372,1880,9550},
		{5825,8517,7769,2405,8204,1060,3603,7025,478,8334,1997,3692,7433,9101,7294,7498,9415,5452,3850,3508,6857,9213,6807,4412,7310,854,5384,686,4978,892,8651,3241,2743,3801,3813,8588,6701,4416,6990,6490,3197,6838,6503,114,8343,5844,8646,8694,65,791,5979,2687,2621,2019,8097,1423,3644,9764,4921,3266,3662,5561,2476,8271,8138,6147,1168,3340,1998,9874,6572,9873,6659,5609,2711,3931,9567,4143,7833,8887},
		{6223,2099,2700,589,4716,8333,1362,5007,2753,2848,4441,8397,7192,8191,4916,9955,6076,3370,6396,6971,3156,248,3911,2488,4930,2458,7183,5455,170,6809,6417,3390,1956,7188,577,7526,2203,968,8164,479,8699,7915,507,6393,4632,1597,7534,3604,618,3280,6061,9793,9238,8347,568,9645,2070,5198,6482,5000,9212,6655,5961,7513,1323,3872,6170,3812,4146,2736,67,3151,5548,2781,9679,7564,5043,8587,1893,4531},
		{5826,3690,6724,2121,9308,6986,8106,6659,2142,1642,7170,2877,5757,6494,8026,6571,8387,9961,6043,9758,9607,6450,8631,8334,7359,5256,8523,2225,7487,1977,9555,8048,5763,2414,4948,4265,2427,8978,8088,8841,9208,9601,5810,9398,8866,9138,4176,5875,7212,3272,6759,5678,7649,4922,5422,1343,8197,3154,3600,687,1028,4579,2084,9467,4492,7262,7296,6538,7657,7134,2077,1505,7332,6890,8964,4879,7603,7400,5973,739},
		{1861,1613,4879,1884,7334,966,2000,7489,2123,4287,1472,3263,4726,9203,1040,4103,6075,6049,330,9253,4062,4268,1635,9960,577,1320,3195,9628,1030,4092,4979,6474,6393,2799,6967,8687,7724,7392,9927,2085,3200,6466,8702,265,7646,8665,7986,7266,4574,6587,612,2724,704,3191,8323,9523,3002,704,5064,3960,8209,2027,2758,8393,4875,4641,9584,6401,7883,7014,768,443,5490,7506,1852,2005,8850,5776,4487,4269},
		{4052,6687,4705,7260,6645,6715,3706,5504,8672,2853,1136,8187,8203,4016,871,1809,1366,4952,9294,5339,6872,2645,6083,7874,3056,5218,7485,8796,7401,3348,2103,426,8572,4163,9171,3176,948,7654,9344,3217,1650,5580,7971,2622,76,2874,880,2034,9929,1546,2659,5811,3754,7096,7436,9694,9960,7415,2164,953,2360,4194,2397,1047,2196,6827,575,784,2675,8821,6802,7972,5996,6699,2134,7577,2887,1412,4349,4380},
		{4629,2234,6240,8132,7592,3181,6389,1214,266,1910,2451,8784,2790,1127,6932,1447,8986,2492,5476,397,889,3027,7641,5083,5776,4022,185,3364,5701,2442,2840,4160,9525,4828,6602,2614,7447,3711,4505,7745,8034,6514,4907,2605,7753,6958,7270,6936,3006,8968,439,2326,4652,3085,3425,9863,5049,5361,8688,297,7580,8777,7916,6687,8683,7141,306,9569,2384,1500,3346,4601,7329,9040,6097,2727,6314,4501,4974,2829},
		{8316,4072,2025,6884,3027,1808,5714,7624,7880,8528,4205,8686,7587,3230,1139,7273,6163,6986,3914,9309,1464,9359,4474,7095,2212,7302,2583,9462,7532,6567,1606,4436,8981,5612,6796,4385,5076,2007,6072,3678,8331,1338,3299,8845,4783,8613,4071,1232,6028,2176,3990,2148,3748,103,9453,538,6745,9110,926,3125,473,5970,8728,7072,9062,1404,1317,5139,9862,6496,6062,3338,464,1600,2532,1088,8232,7739,8274,3873},
		{2341,523,7096,8397,8301,6541,9844,244,4993,2280,7689,4025,4196,5522,7904,6048,2623,9258,2149,9461,6448,8087,7245,1917,8340,7127,8466,5725,6996,3421,5313,512,9164,9837,9794,8369,4185,1488,7210,1524,1016,4620,9435,2478,7765,8035,697,6677,3724,6988,5853,7662,3895,9593,1185,4727,6025,5734,7665,3070,138,8469,6748,6459,561,7935,8646,2378,462,7755,3115,9690,8877,3946,2728,8793,244,6323,8666,4271},
		{6430,2406,8994,56,1267,3826,9443,7079,7579,5232,6691,3435,6718,5698,4144,7028,592,2627,217,734,6194,8156,9118,58,2640,8069,4127,3285,694,3197,3377,4143,4802,3324,8134,6953,7625,3598,3584,4289,7065,3434,2106,7132,5802,7920,9060,7531,3321,1725,1067,3751,444,5503,6785,7937,6365,4803,198,6266,8177,1470,6390,1606,2904,7555,9834,8667,2033,1723,5167,1666,8546,8152,473,4475,6451,7947,3062,3281},
		{2810,3042,7759,1741,2275,2609,7676,8640,4117,1958,7500,8048,1757,3954,9270,1971,4796,2912,660,5511,3553,1012,5757,4525,6084,7198,8352,5775,7726,8591,7710,9589,3122,4392,6856,5016,749,2285,3356,7482,9956,7348,2599,8944,495,3462,3578,551,4543,7207,7169,7796,1247,4278,6916,8176,3742,8385,2310,1345,8692,2667,4568,1770,8319,3585,4920,3890,4928,7343,5385,9772,7947,8786,2056,9266,3454,2807,877,2660},
		{6206,8252,5928,5837,4177,4333,207,7934,5581,9526,8906,1498,8411,2984,5198,5134,2464,8435,8514,8674,3876,599,5327,826,2152,4084,2433,9327,9697,4800,2728,3608,3849,3861,3498,9943,1407,3991,7191,9110,5666,8434,4704,6545,5944,2357,1163,4995,9619,6754,4200,9682,6654,4862,4744,5953,6632,1054,293,9439,8286,2255,696,8709,1533,1844,6441,430,1999,6063,9431,7018,8057,2920,6266,6799,356,3597,4024,6665},
		{3847,6356,8541,7225,2325,2946,5199,469,5450,7508,2197,9915,8284,7983,6341,3276,3321,16,1321,7608,5015,3362,8491,6968,6818,797,156,2575,706,9516,5344,5457,9210,5051,8099,1617,9951,7663,8253,9683,2670,1261,4710,1068,8753,4799,1228,2621,3275,6188,4699,1791,9518,8701,5932,4275,6011,9877,2933,4182,6059,2930,6687,6682,9771,654,9437,3169,8596,1827,5471,8909,2352,123,4394,3208,8756,5513,6917,2056},
		{5458,8173,3138,3290,4570,4892,3317,4251,9699,7973,1163,1935,5477,6648,9614,5655,9592,975,9118,2194,7322,8248,8413,3462,8560,1907,7810,6650,7355,2939,4973,6894,3933,3784,3200,2419,9234,4747,2208,2207,1945,2899,1407,6145,8023,3484,5688,7686,2737,3828,3704,9004,5190,9740,8643,8650,5358,4426,1522,1707,3613,9887,6956,2447,2762,833,1449,9489,2573,1080,4167,3456,6809,2466,227,7125,2759,6250,6472,8089},
		{3266,7025,9756,3914,1265,9116,7723,9788,6805,5493,2092,8688,6592,9173,4431,4028,6007,7131,4446,4815,3648,6701,759,3312,8355,4485,4187,5188,8746,7759,3528,2177,5243,8379,3838,7233,4607,9187,7216,2190,6967,2920,6082,7910,5354,3609,8958,6949,7731,494,8753,8707,1523,4426,3543,7085,647,6771,9847,646,5049,824,8417,5260,2730,5702,2513,9275,4279,2767,8684,1165,9903,4518,55,9682,8963,6005,2102,6523},
		{1998,8731,936,1479,5259,7064,4085,91,7745,7136,3773,3810,730,8255,2705,2653,9790,6807,2342,355,9344,2668,3690,2028,9679,8102,574,4318,6481,9175,5423,8062,2867,9657,7553,3442,3920,7430,3945,7639,3714,3392,2525,4995,4850,2867,7951,9667,486,9506,9888,781,8866,1702,3795,90,356,1483,4200,2131,6969,5931,486,6880,4404,1084,5169,4910,6567,8335,4686,5043,2614,3352,2667,4513,6472,7471,5720,1616},
		{8878,1613,1716,868,1906,2681,564,665,5995,2474,7496,3432,9491,9087,8850,8287,669,823,347,6194,2264,2592,7871,7616,8508,4827,760,2676,4660,4881,7572,3811,9032,939,4384,929,7525,8419,5556,9063,662,8887,7026,8534,3111,1454,2082,7598,5726,6687,9647,7608,73,3014,5063,670,5461,5631,3367,9796,8475,7908,5073,1565,5008,5295,4457,1274,4788,1728,338,600,8415,8535,9351,7750,6887,5845,1741,125},
		{3637,6489,9634,9464,9055,2413,7824,9517,7532,3577,7050,6186,6980,9365,9782,191,870,2497,8498,2218,2757,5420,6468,586,3320,9230,1034,1393,9886,5072,9391,1178,8464,8042,6869,2075,8275,3601,7715,9470,8786,6475,8373,2159,9237,2066,3264,5000,679,355,3069,4073,494,2308,5512,4334,9438,8786,8637,9774,1169,1949,6594,6072,4270,9158,7916,5752,6794,9391,6301,5842,3285,2141,3898,8027,4310,8821,7079,1307},
		{8497,6681,4732,7151,7060,5204,9030,7157,833,5014,8723,3207,9796,9286,4913,119,5118,7650,9335,809,3675,2597,5144,3945,5090,8384,187,4102,1260,2445,2792,4422,8389,9290,50,1765,1521,6921,8586,4368,1565,5727,7855,2003,4834,9897,5911,8630,5070,1330,7692,7557,7980,6028,5805,9090,8265,3019,3802,698,9149,5748,1965,9658,4417,5994,5584,8226,2937,272,5743,1278,5698,8736,2595,6475,5342,6596,1149,6920},
		{8188,8009,9546,6310,8772,2500,9846,6592,6872,3857,1307,8125,7042,1544,6159,2330,643,4604,7899,6848,371,8067,2062,3200,7295,1857,9505,6936,384,2193,2190,301,8535,5503,1462,7380,5114,4824,8833,1763,4974,8711,9262,6698,3999,2645,6937,7747,1128,2933,3556,7943,2885,3122,9105,5447,418,2899,5148,3699,9021,9501,597,4084,175,1621,1,1079,6067,5812,4326,9914,6633,5394,4233,6728,9084,1864,5863,1225},
		{9935,8793,9117,1825,9542,8246,8437,3331,9128,9675,6086,7075,319,1334,7932,3583,7167,4178,1726,7720,695,8277,7887,6359,5912,1719,2780,8529,1359,2013,4498,8072,1129,9998,1147,8804,9405,6255,1619,2165,7491,1,8882,7378,3337,503,5758,4109,3577,985,3200,7615,8058,5032,1080,6410,6873,5496,1466,2412,9885,5904,4406,3605,8770,4361,6205,9193,1537,9959,214,7260,9566,1685,100,4920,7138,9819,5637,976},
		{3466,9854,985,1078,7222,8888,5466,5379,3578,4540,6853,8690,3728,6351,7147,3134,6921,9692,857,3307,4998,2172,5783,3931,9417,2541,6299,13,787,2099,9131,9494,896,8600,1643,8419,7248,2660,2609,8579,91,6663,5506,7675,1947,6165,4286,1972,9645,3805,1663,1456,8853,5705,9889,7489,1107,383,4044,2969,3343,152,7805,4980,9929,5033,1737,9953,7197,9158,4071,1324,473,9676,3984,9680,3606,8160,7384,5432},
		{1005,4512,5186,3953,2164,3372,4097,3247,8697,3022,9896,4101,3871,6791,3219,2742,4630,6967,7829,5991,6134,1197,1414,8923,8787,1394,8852,5019,7768,5147,8004,8825,5062,9625,7988,1110,3992,7984,9966,6516,6251,8270,421,3723,1432,4830,6935,8095,9059,2214,6483,6846,3120,1587,6201,6691,9096,9627,6671,4002,3495,9939,7708,7465,5879,6959,6634,3241,3401,2355,9061,2611,7830,3941,2177,2146,5089,7079,519,6351},
		{7280,8586,4261,2831,7217,3141,9994,9940,5462,2189,4005,6942,9848,5350,8060,6665,7519,4324,7684,657,9453,9296,2944,6843,7499,7847,1728,9681,3906,6353,5529,2822,3355,3897,7724,4257,7489,8672,4356,3983,1948,6892,7415,4153,5893,4190,621,1736,4045,9532,7701,3671,1211,1622,3176,4524,9317,7800,5638,6644,6943,5463,3531,2821,1347,5958,3436,1438,2999,994,850,4131,2616,1549,3465,5946,690,9273,6954,7991},
		{9517,399,3249,2596,7736,2142,1322,968,7350,1614,468,3346,3265,7222,6086,1661,5317,2582,7959,4685,2807,2917,1037,5698,1529,3972,8716,2634,3301,3412,8621,743,8001,4734,888,7744,8092,3671,8941,1487,5658,7099,2781,99,1932,4443,4756,4652,9328,1581,7855,4312,5976,7255,6480,3996,2748,1973,9731,4530,2790,9417,7186,5303,3557,351,7182,9428,1342,9020,7599,1392,8304,2070,9138,7215,2008,9937,1106,7110},
		{7444,769,9688,632,1571,6820,8743,4338,337,3366,3073,1946,8219,104,4210,6986,249,5061,8693,7960,6546,1004,8857,5997,9352,4338,6105,5008,2556,6518,6694,4345,3727,7956,20,3954,8652,4424,9387,2035,8358,5962,5304,5194,8650,8282,1256,1103,2138,6679,1985,3653,2770,2433,4278,615,2863,1715,242,3790,2636,6998,3088,1671,2239,957,5411,4595,6282,2881,9974,2401,875,7574,2987,4587,3147,6766,9885,2965},
		{3287,3016,3619,6818,9073,6120,5423,557,2900,2015,8111,3873,1314,4189,1846,4399,7041,7583,2427,2864,3525,5002,2069,748,1948,6015,2684,438,770,8367,1663,7887,7759,1885,157,7770,4520,4878,3857,1137,3525,3050,6276,5569,7649,904,4533,7843,2199,5648,7628,9075,9441,3600,7231,2388,5640,9096,958,3058,584,5899,8150,1181,9616,1098,8162,6819,8171,1519,1140,7665,8801,2632,1299,9192,707,9955,2710,7314},
		{1772,2963,7578,3541,3095,1488,7026,2634,6015,4633,4370,2762,1650,2174,909,8158,2922,8467,4198,4280,9092,8856,8835,5457,2790,8574,9742,5054,9547,4156,7940,8126,9824,7340,8840,6574,3547,1477,3014,6798,7134,435,9484,9859,3031,4,1502,4133,1738,1807,4825,463,6343,9701,8506,9822,9555,8688,8168,3467,3234,6318,1787,5591,419,6593,7974,8486,9861,6381,6758,194,3061,4315,2863,4665,3789,2201,1492,4416},
		{126,8927,6608,5682,8986,6867,1715,6076,3159,788,3140,4744,830,9253,5812,5021,7616,8534,1546,9590,1101,9012,9821,8132,7857,4086,1069,7491,2988,1579,2442,4321,2149,7642,6108,250,6086,3167,24,9528,7663,2685,1220,9196,1397,5776,1577,1730,5481,977,6115,199,6326,2183,3767,5928,5586,7561,663,8649,9688,949,5913,9160,1870,5764,9887,4477,6703,1413,4995,5494,7131,2192,8969,7138,3997,8697,646,1028},
		{8074,1731,8245,624,4601,8706,155,8891,309,2552,8208,8452,2954,3124,3469,4246,3352,1105,4509,8677,9901,4416,8191,9283,5625,7120,2952,8881,7693,830,4580,8228,9459,8611,4499,1179,4988,1394,550,2336,6089,6872,269,7213,1848,917,6672,4890,656,1478,6536,3165,4743,4990,1176,6211,7207,5284,9730,4738,1549,4986,4942,8645,3698,9429,1439,2175,6549,3058,6513,1574,6988,8333,3406,5245,5431,7140,7085,6407},
		{7845,4694,2530,8249,290,5948,5509,1588,5940,4495,5866,5021,4626,3979,3296,7589,4854,1998,5627,3926,8346,6512,9608,1918,7070,4747,4182,2858,2766,4606,6269,4107,8982,8568,9053,4244,5604,102,2756,727,5887,2566,7922,44,5986,621,1202,374,6988,4130,3627,6744,9443,4568,1398,8679,397,3928,9159,367,2917,6127,5788,3304,8129,911,2669,1463,9749,264,4478,8940,1109,7309,2462,117,4692,7724,225,2312},
		{4164,3637,2000,941,8903,39,3443,7172,1031,3687,4901,8082,4945,4515,7204,9310,9349,9535,9940,218,1788,9245,2237,1541,5670,6538,6047,5553,9807,8101,1925,8714,445,8332,7309,6830,5786,5736,7306,2710,3034,1838,7969,6318,7912,2584,2080,7437,6705,2254,7428,820,782,9861,7596,3842,3631,8063,5240,6666,394,4565,7865,4895,9890,6028,6117,4724,9156,4473,4552,602,470,6191,4927,5387,884,3146,1978,3000},
		{4258,6880,1696,3582,5793,4923,2119,1155,9056,9698,6603,3768,5514,9927,9609,6166,6566,4536,4985,4934,8076,9062,6741,6163,7399,4562,2337,5600,2919,9012,8459,1308,6072,1225,9306,8818,5886,7243,7365,8792,6007,9256,6699,7171,4230,7002,8720,7839,4533,1671,478,7774,1607,2317,5437,4705,7886,4760,6760,7271,3081,2997,3088,7675,6208,3101,6821,6840,122,9633,4900,2067,8546,4549,2091,7188,5605,8599,6758,5229},
		{7854,5243,9155,3556,8812,7047,2202,1541,5993,4600,4760,713,434,7911,7426,7414,8729,322,803,7960,7563,4908,6285,6291,736,3389,9339,4132,8701,7534,5287,3646,592,3065,7582,2592,8755,6068,8597,1982,5782,1894,2900,6236,4039,6569,3037,5837,7698,700,7815,2491,7272,5878,3083,6778,6639,3589,5010,8313,2581,6617,5869,8402,6808,2951,2321,5195,497,2190,6187,1342,1316,4453,7740,4154,2959,1781,1482,8256},
		{7178,2046,4419,744,8312,5356,6855,8839,319,2962,5662,47,6307,8662,68,4813,567,2712,9931,1678,3101,8227,6533,4933,6656,92,5846,4780,6256,6361,4323,9985,1231,2175,7178,3034,9744,6155,9165,7787,5836,9318,7860,9644,8941,6480,9443,8188,5928,161,6979,2352,5628,6991,1198,8067,5867,6620,3778,8426,2994,3122,3124,6335,3918,8897,2655,9670,634,1088,1576,8935,7255,474,8166,7417,9547,2886,5560,3842},
		{6957,3111,26,7530,7143,1295,1744,6057,3009,1854,8098,5405,2234,4874,9447,2620,9303,27,7410,969,40,2966,5648,7596,8637,4238,3143,3679,7187,690,9980,7085,7714,9373,5632,7526,6707,3951,9734,4216,2146,3602,5371,6029,3039,4433,4855,4151,1449,3376,8009,7240,7027,4602,2947,9081,4045,8424,9352,8742,923,2705,4266,3232,2264,6761,363,2651,3383,7770,6730,7856,7340,9679,2158,610,4471,4608,910,6241},
		{4417,6756,1013,8797,658,8809,5032,8703,7541,846,3357,2920,9817,1745,9980,7593,4667,3087,779,3218,6233,5568,4296,2289,2654,7898,5021,9461,5593,8214,9173,4203,2271,7980,2983,5952,9992,8399,3468,1776,3188,9314,1720,6523,2933,621,8685,5483,8986,6163,3444,9539,4320,155,3992,2828,2150,6071,524,2895,5468,8063,1210,3348,9071,4862,483,9017,4097,6186,9815,3610,5048,1644,1003,9865,9332,2145,1944,2213},
		{9284,3803,4920,1927,6706,4344,7383,4786,9890,2010,5228,1224,3158,6967,8580,8990,8883,5213,76,8306,2031,4980,5639,9519,7184,5645,7769,3259,8077,9130,1317,3096,9624,3818,1770,695,2454,947,6029,3474,9938,3527,5696,4760,7724,7738,2848,6442,5767,6845,8323,4131,2859,7595,2500,4815,3660,9130,8580,7016,8231,4391,8369,3444,4069,4021,556,6154,627,2778,1496,4206,6356,8434,8491,3816,8231,3190,5575,1015},
		{3787,7572,1788,6803,5641,6844,1961,4811,8535,9914,9999,1450,8857,738,4662,8569,6679,2225,7839,8618,286,2648,5342,2294,3205,4546,176,8705,3741,6134,8324,8021,7004,5205,7032,6637,9442,5539,5584,4819,5874,5807,8589,6871,9016,983,1758,3786,1519,6241,185,8398,495,3370,9133,3051,4549,9674,7311,9738,3316,9383,2658,2776,9481,7558,619,3943,3324,6491,4933,153,9738,4623,912,3595,7771,7939,1219,4405},
		{2650,3883,4154,5809,315,7756,4430,1788,4451,1631,6461,7230,6017,5751,138,588,5282,2442,9110,9035,6349,2515,1570,6122,4192,4174,3530,1933,4186,4420,4609,5739,4135,2963,6308,1161,8809,8619,2796,3819,6971,8228,4188,1492,909,8048,2328,6772,8467,7671,9068,2226,7579,6422,7056,8042,3296,2272,3006,2196,7320,3238,3490,3102,37,1293,3212,4767,5041,8773,5794,4456,6174,7279,7054,2835,7053,9088,790,6640},
		{3101,1057,7057,3826,6077,1025,2955,1224,1114,6729,5902,4698,6239,7203,9423,1804,4417,6686,1426,6941,8071,1029,4985,9010,6122,6597,1622,1574,3513,1684,7086,5505,3244,411,9638,4150,907,9135,829,981,1707,5359,8781,9751,5,9131,3973,7159,1340,6955,7514,7993,6964,8198,1933,2797,877,3993,4453,8020,9349,8646,2779,8679,2961,3547,3374,3510,1129,3568,2241,2625,9138,5974,8206,7669,7678,1833,8700,4480},
		{4865,9912,8038,8238,782,3095,8199,1127,4501,7280,2112,2487,3626,2790,9432,1475,6312,8277,4827,2218,5806,7132,8752,1468,7471,6386,739,8762,8323,8120,5169,9078,9058,3370,9560,7987,8585,8531,5347,9312,1058,4271,1159,5286,5404,6925,8606,9204,7361,2415,560,586,4002,2644,1927,2824,768,4409,2942,3345,1002,808,4941,6267,7979,5140,8643,7553,9438,7320,4938,2666,4609,2778,8158,6730,3748,3867,1866,7181},
		{171,3771,7134,8927,4778,2913,3326,2004,3089,7853,1378,1729,4777,2706,9578,1360,5693,3036,1851,7248,2403,2273,8536,6501,9216,613,9671,7131,7719,6425,773,717,8803,160,1114,7554,7197,753,4513,4322,8499,4533,2609,4226,8710,6627,644,9666,6260,4870,5744,7385,6542,6203,7703,6130,8944,5589,2262,6803,6381,7414,6888,5123,7320,9392,9061,6780,322,8975,7050,5089,1061,2260,3199,1150,1865,5386,9699,6501},
		{3744,8454,6885,8277,919,1923,4001,6864,7854,5519,2491,6057,8794,9645,1776,5714,9786,9281,7538,6916,3215,395,2501,9618,4835,8846,9708,2813,3303,1794,8309,7176,2206,1602,1838,236,4593,2245,8993,4017,10,8215,6921,5206,4023,5932,6997,7801,262,7640,3107,8275,4938,7822,2425,3223,3886,2105,8700,9526,2088,8662,8034,7004,5710,2124,7164,3574,6630,9980,4242,2901,9471,1491,2117,4562,1130,9086,4117,6698},
		{2810,2280,2331,1170,4554,4071,8387,1215,2274,9848,6738,1604,7281,8805,439,1298,8318,7834,9426,8603,6092,7944,1309,8828,303,3157,4638,4439,9175,1921,4695,7716,1494,1015,1772,5913,1127,1952,1950,8905,4064,9890,385,9357,7945,5035,7082,5369,4093,6546,5187,5637,2041,8946,1758,7111,6566,1027,1049,5148,7224,7248,296,6169,375,1656,7993,2816,3717,4279,4675,1609,3317,42,6201,3100,3144,163,9530,4531},
		{7096,6070,1009,4988,3538,5801,7149,3063,2324,2912,7911,7002,4338,7880,2481,7368,3516,2016,7556,2193,1388,3865,8125,4637,4096,8114,750,3144,1938,7002,9343,4095,1392,4220,3455,6969,9647,1321,9048,1996,1640,6626,1788,314,9578,6630,2813,6626,4981,9908,7024,4355,3201,3521,3864,3303,464,1923,595,9801,3391,8366,8084,9374,1041,8807,9085,1892,9431,8317,9016,9221,8574,9981,9240,5395,2009,6310,2854,9255},
		{8830,3145,2960,9615,8220,6061,3452,2918,6481,9278,2297,3385,6565,7066,7316,5682,107,7646,4466,68,1952,9603,8615,54,7191,791,6833,2560,693,9733,4168,570,9127,9537,1925,8287,5508,4297,8452,8795,6213,7994,2420,4208,524,5915,8602,8330,2651,8547,6156,1812,6271,7991,9407,9804,1553,6866,1128,2119,4691,9711,8315,5879,9935,6900,482,682,4126,1041,428,6247,3720,5882,7526,2582,4327,7725,3503,2631},
		{2738,9323,721,7434,1453,6294,2957,3786,5722,6019,8685,4386,3066,9057,6860,499,5315,3045,5194,7111,3137,9104,941,586,3066,755,4177,8819,7040,5309,3583,3897,4428,7788,4721,7249,6559,7324,825,7311,3760,6064,6070,9672,4882,584,1365,9739,9331,5783,2624,7889,1604,1303,1555,7125,8312,425,8936,3233,7724,1480,403,7440,1784,1754,4721,1569,652,3893,4574,5692,9730,4813,9844,8291,9199,7101,3391,8914},
		{6044,2928,9332,3328,8588,447,3830,1176,3523,2705,8365,6136,5442,9049,5526,8575,8869,9031,7280,706,2794,8814,5767,4241,7696,78,6570,556,5083,1426,4502,3336,9518,2292,1885,3740,3153,9348,9331,8051,2759,5407,9028,7840,9255,831,515,2612,9747,7435,8964,4971,2048,4900,5967,8271,1719,9670,2810,6777,1594,6367,6259,8316,3815,1689,6840,9437,4361,822,9619,3065,83,6344,7486,8657,8228,9635,6932,4864},
		{8478,4777,6334,4678,7476,4963,6735,3096,5860,1405,5127,7269,7793,4738,227,9168,2996,8928,765,733,1276,7677,6258,1528,9558,3329,302,8901,1422,8277,6340,645,9125,8869,5952,141,8141,1816,9635,4025,4184,3093,83,2344,2747,9352,7966,1206,1126,1826,218,7939,2957,2729,810,8752,5247,4174,4038,8884,7899,9567,301,5265,5752,7524,4381,1669,3106,8270,6228,6373,754,2547,4240,2313,5514,3022,1040,9738},
		{2265,8192,1763,1369,8469,8789,4836,52,1212,6690,5257,8918,6723,6319,378,4039,2421,8555,8184,9577,1432,7139,8078,5452,9628,7579,4161,7490,5159,8559,1011,81,478,5840,1964,1334,6875,8670,9900,739,1514,8692,522,9316,6955,1345,8132,2277,3193,9773,3923,4177,2183,1236,6747,6575,4874,6003,6409,8187,745,8776,9440,7543,9825,2582,7381,8147,7236,5185,7564,6125,218,7991,6394,391,7659,7456,5128,5294},
		{2132,8992,8160,5782,4420,3371,3798,5054,552,5631,7546,4716,1332,6486,7892,7441,4370,6231,4579,2121,8615,1145,9391,1524,1385,2400,9437,2454,7896,7467,2928,8400,3299,4025,7458,4703,7206,6358,792,6200,725,4275,4136,7390,5984,4502,7929,5085,8176,4600,119,3568,76,9363,6943,2248,9077,9731,6213,5817,6729,4190,3092,6910,759,2682,8380,1254,9604,3011,9291,5329,9453,9746,2739,6522,3765,5634,1113,5789},
		{5304,5499,564,2801,679,2653,1783,3608,7359,7797,3284,796,3222,437,7185,6135,8571,2778,7488,5746,678,6140,861,7750,803,9859,9918,2425,3734,2698,9005,4864,9818,6743,2475,132,9486,3825,5472,919,292,4411,7213,7699,6435,9019,6769,1388,802,2124,1345,8493,9487,8558,7061,8777,8833,2427,2238,5409,4957,8503,3171,7622,5779,6145,2417,5873,5563,5693,9574,9491,1937,7384,4563,6842,5432,2751,3406,7981},
	};
};

std::string Euler082()
{
    //260324
    p082 c;
    return c.run();
}

void visit_sideway(long long N, long long i, long long j, mat::matrix& mminpaths, mat::matrix& mnumber, mat::matrix& mvisit)
{
    long long t1 = 0;
    long long t2 = 0;
    if (mvisit(i,j) == 0)
    {
        if (mnumber(i,j) == 0) {mvisit(i,j) = 1; return;}

        mvisit(i,j) = 1;
        if (i+1<=N-1)                if (mvisit(i+1,j  ) == 0) visit_sideway(N, i+1, j,   mminpaths, mnumber, mvisit);
        if ((i+0<=N-1) && (j+1<=N))  if (mvisit(i+0,j+1) == 0) visit_sideway(N, i+0, j+1, mminpaths, mnumber, mvisit);

        t1 = 0; if (i+1<=N-1)                  t1 = mminpaths(i+1,j );
        t2 = 0; if ((i+0<=N-1) && (j+1<=N-1))  t2 = mminpaths(i+0,j+1);
        if (t1==0) mminpaths(i, j) = mnumber(i, j)  + t2;
        else if (t2==0) mminpaths(i, j) = mnumber(i, j)  + t1;
        else mminpaths(i, j) = mnumber(i, j)  + std::min(t1,t2);
    }
}


long long Euler081(long long N)
{
    //427337
    int NROWS = (int)N;
    mat::matrix mnumber(NROWS, NROWS, 0);
    mat::matrix mminpaths(NROWS, NROWS, 0);
    mat::matrix mvisit(NROWS, NROWS, 0);
    std::string s;

    #ifdef _WIN32
    std::ifstream is("p081_matrix.txt", std::ifstream::in);
    #else
    std::ifstream is("../Euler/p081_matrix.txt", std::ifstream::in);
    #endif
    if (is)
    {
        // get length of file:
        is.seekg(0, is.end);
        int length = (int)is.tellg();
        is.seekg(0, is.beg);

        char* buffer = new char[length];
        //std::cout << "Reading " << length << " characters... " << std::endl;;
        is.read(buffer, length);
        if (is)
        {
            int pos = 0;
            int i = 0;
            int j = 0;

            while(true)
            {
                {
                    if ((buffer[pos]==',') || (buffer[pos]=='\n'))
                    {
                        if (s.size()>0)
                        {
                            mnumber(i, j) = to_long(s);
                            //if (i<3)std::cout << i << " " << j << " " << mnumber(i, j) << std::endl;
                            s.clear();
                            j++;
                            if (j>=NROWS) {i++;j=0;}
                        }
                    }
                    else if (buffer[pos]!=',') s=s+buffer[pos];
                    pos += 1;
                }
                if (pos >= length) break;
            }
            if (s.size()>0)
                mnumber(NROWS-1, NROWS-1) = to_long(s);
        }
        is.close();
    }
    visit_sideway(NROWS, 0, 0, mminpaths, mnumber, mvisit);
    return mminpaths(0, 0);
}

long long Euler080(long long N)
{
    //40886 including digit before decimal!
    dec101_t s;
    long long sint;
    dec101_t sum = 0;
    std::string str;
    int a;

    for (long long i = 2; i <= N; i++)
    {
        s = sqrt(dec101_t(i));
        sint = (long long)s;
        if (sint * sint != i )
        {
            std::stringstream os;
            os << std::setprecision(1000) << s;
            str = os.str();

            int cnt = 0;
            for (size_t j = 0; j < str.size(); j++)
            {
                if (str.at(j) != '.')
                if (cnt < 100)
                {
                    a = str.at(j) - '0';
                    sum += a;
                    cnt++;
                }
            }
        }
    }
    return (long long)sum;
}


long long Euler079(long long N)
{
    // 73162890
    std::vector<bool> vbefore[10];
    for (long long i = 0; i < 10; i++)
    {
        for (long long j = 0; j < 10; j++)
        {
            vbefore[i].push_back(false);
        }
    }

    // READ...
    int num[50 * 3];
    std::ifstream is("p079_keylog.txt", std::ifstream::in);
    if (is)
    {
        // get length of file:
        is.seekg(0, is.end);
        int length = (int)is.tellg();
        is.seekg(0, is.beg);

        char* buffer = new char[length];

        //std::cout << "Reading " << length << " characters... " << std::endl;;
        // read data as a block:
        is.read(buffer, length);
        if (is)
        {
            int pos = 0;
            for (int i = 0; i < 50; i++)
            {
                for (int j = 0; j < 4; j++)
                {
                    if (j!= 3)
                        num[3 * i + j] = (int)(buffer[pos] - '0');
                    pos += 1;
                    if (j != 3)
                        std::cout << i << " " << j << " " << num[3 * i + j] << std::endl;
                }
            }
        }
        else
        {
            std::cout << "NO FILE" << std::endl;
            return 0;
        }
        is.close();
    }

    //319
    int a; int b; int c;
    std::map<int, bool> vmap;
    for (int i = 0; i < 50; i++)
    {
        // 319 3 before 1 1 before 9 3x1x9
        {
            a = num[3 * i];
            b = num[3 * i + 1];
            c = num[3 * i + 2];
            vbefore[a][b] = true;
            vbefore[a][c] = true;
            vbefore[b][c] = true;

            vmap[a] = true;
            vmap[b] = true;
            vmap[c] = true;
        }
    }

    int posj;
    bool ok;
    for (long long i = 10000000; i < N; i++) // file has minimum 8 digits
    {
        ok = true;
        auto v = digits10(i, true);
        if (v.size() >= 8)
        {
            // assume not repeating key...
            if (all_digit_unique(v))
            {
                for (long long  ii = 0; ii < (long long)v.size(); ii++)
                {
                    if (ok)
                    {
                        for (long long j = 0; j < 10; j++)
                        {
                            if ((vmap[(int)j] == false) && (count_digit((int)j, v) > 0))
                            {
                                ok = false;
                                break;
                            }

                            if (vbefore[v[ii]][j] == true)
                            {
                                posj = -1;
                                for (long long d = 0; d < (long long)v.size(); d++)
                                {
                                    if (v[d] == j)
                                    {
                                        posj = (int)d;
                                        break;
                                    }
                                }
                                if (posj == -1)
                                {
                                    ok = false;
                                    break;
                                }
                                if (posj <= ii)
                                {
                                    ok = false;
                                    break;
                                }
                            }
                        }
                    }
                }
                if (ok) return i;
            }
        }
    }

    // First is the one before all other recurselively ...
    return 0;
}

//
//The multiplicative inverse of its generating function is the Euler function;
//by Euler's pentagonal number theorem this function is an alternating sum of pentagonal number powers of its argument.
//p(n)=p(n-1)+p(n-2)-p(n-5)-p(n-7)...
long long pentagonal_number(long  long k)
{
    long long r = (3 * k * k - k) / 2;
    return r;
}

long long Euler078(long long N)
{
    //55374
    std::map<long long, uinteger_t> partionmap;
    std::vector<long long> vpentagonal(N+1, 0);
    long long cnt = 0;
    for (long long i = 1; i < N/2; i++)
    {
        vpentagonal[cnt] = pentagonal_number(i); // 1,2   5,7   12,15   22.. 35.. ..
        cnt++;
        vpentagonal[cnt] = pentagonal_number(-1 * i);
        cnt++;
    }

    partionmap[0] = 1;
    partionmap[1] = 1;
    partionmap[2] = 2;
    partionmap[3] = 3;
    partionmap[4] = 5;
    partionmap[5] = 7; // 11 15 22 30 42 56
    long long pn;
    uinteger_t sum;

    for (long long i = 3; i <= N; i++)
    {
        sum = 0;
        int plus = 1;
        int n = 0;
        for (long long p = 0; p < (long long)vpentagonal.size(); p++)
        {
            if (i - vpentagonal[p] >= 0)
            {
                pn = i - vpentagonal[p];
                if (plus == 1)
                    sum += partionmap[pn];
                else
                    sum -= partionmap[pn];
                n++;
                if (n % 2 == 0) plus = 1 - plus;
            }
            else
            {
                break;
            }
        }
        std::cout << i << " " << sum << " " << std::endl;

        partionmap[i] = sum;
        if (sum % 1000000 == 0) return i;
    }
    return 0;
}

long long Euler077(long long N)
{
    // 71
    // Recursive count same idea of Euler076
    long long r = 0;
    constexpr int MAXPRIME = 500;; // guess
    std::vector<long long> primes = list_prime(MAXPRIME);
    long long ways[MAXPRIME] = { 0 };
    ways[0] = 1;
    for (long long i = 0; i < (long long)primes.size(); i++)
    {
        for (long long j = 0; j + primes[i] < MAXPRIME; j++)
        {
            ways[j + primes[i]] += ways[j];
        }
        for (long long j = (i > 0 ? primes[i - 1] + 1 : 0); j <= primes[i]; j++)
        {
            if (ways[j] >= N)
            {
                r = j;
                i = primes.size();
                break;
            }
        }
    }
    return r;
}

long long Euler076(long long N)
{
    //190569291
    //EX: P(7)= P(6) + P(5) +... + P(1)
    int* count = new int[N+1];
    for (int m = 0; m < N+1; m++) count[m] = 0;

    count[0] = 1;
    for (int i = 1; i < N; i++) {
        for (int j = 0; j + i < N + 1; j++) {
            count[i + j] += count[j];
        }
    }
    auto r = count[N];;
    delete []count;
    return r;
}

long long Euler075(long long N)
{
    //161667
    //Need Euclid's formula: x�+y�=(2mn)�+(m�-n�)�=(m�+n�)�=z�. L=2mn+m2-n2+m2+n2 = 2 * m * (m + n)
    int MAXL = (int)(N+1);
    char* list = new char[MAXL];
    for (int m = 0; m < MAXL; m++) list[m] = 0;

    int crossbound = (int)sqrt((double)MAXL);
    for (int m = 2; m <= crossbound; m++)
    {
        std::cout << m << std::endl;
        for (int n = 1; n < m; n += 2)
        {
            //(m-n) must be odd
            if ((m - n) % 2 == 0)
            {
                n--;
                continue;
            }
            int l = 2 * m * (m + n); //Reduced form of Euclid's formula
            if (l >= MAXL) break;
            if (list[l] <= 1 && gcd(m, n) == 1)
            {
                for (int i = l; i < MAXL; i += l)
                    list[i]++;
            }
        }
    }
    int res = 0;
    for (int i = 1; i < MAXL; i++)
        if (list[i] == 1)
            res++;

    delete[]list;
    return res;
}


long long Euler074(long long N)
{
    //402
    long long cnt = 0;
    long long sum = 0;
    bool ok;
    std::vector<int> va;
    std::vector<int> vsum;
    std::map<std::vector<int>, bool> vmap;

    for (long long a = 1; a < N; a++)
    {
        vmap.clear();
        va = digits10(a);
        vmap[va] = true;

        ok = true;
        while (ok)
        {
            sum = 0;
            for (long long b = 0; b < (long long)va.size(); b++)
            {
                sum += longfact(va[b]);
            }
            vsum = digits10(sum);
            if (vmap.find(vsum) == vmap.end())
            {
                vmap[vsum] = true;
            }
            else
            {
                if (vmap.size() == 60)
                {
                    cnt++;
                }
                ok = false;
                break;
            }
            std::swap(va, vsum);
            vsum.clear();
        }
        if (a % 1000 == 0) std::cout << a << " " << cnt << " " << std::endl;
    }
    return cnt;
}


long long Farey(long long N)
{
    long long cnt = 1;
    for(long long a = 1; a<= N;a++)
    {
        cnt += (long long)phi(a);
    }
    return cnt;
}

long long HCF(long long  a, long long  b)
{
    long long  t;
    while (b != 0)
    {
        t = a;
        a = b;
        b = t % b;
    }
    return a;
}

long long Euler073(long long N)
{
    //Farey sequence
    //7295372
    long long cnt = 0;
    long long low; long long high;

    for(long long a = 4; a<= N;a++)
    {
        low = a/3; low++;
        if (a%3 == 0) low--;
        high = a/2;
        if (a%2 == 0) high--;
        if (high >= low)
        {
            for(long long b = low; b<= high; b++)
            {
                // no common divisor...
                if (HCF(a, b) == 1)
                {
                    cnt++;
                }
            }
        }
    }
    return cnt;
}

long long Euler072(long long N)
{
    //303963552391
    long long cnt = 0;
    for(long long a = 2; a< N;a++)
    {
        if (is_prime(a))
        {
            // 1/5 2/5 ...
            cnt += (a-1);
        }
        else
        {
            cnt += 1;
            cnt += (long long)phi(a);
        }
        if (a%1000 == 0) std::cout << a << " " << cnt << " " << std::endl;
    }
    return cnt;
}

long long Euler071(long long N)
{
    //428570
    RationalNumber r;
    RationalNumber r37 = RationalNumber(3, 7);
    RationalNumber rNearestLess =  RationalNumber(2, 5);
    for(long long a = 1; a< N;a++)
    {
        //b/a < 3/7
        long long bs = (long long) ( - 2 + ((double)a) * 3.0 / 7.0 );
        for(long long b = bs; b < a;b++)
        {
            r = RationalNumber(b, a);
            if ((r < r37) && (r > rNearestLess))
            {
                rNearestLess = r;
                std::cout << a << " " << b << " " << std::endl;
            }
            if (r > r37) break;
        }
    }
    return rNearestLess.getM().toInt(rNearestLess.getM().getNumber());
    //return (long long)rNearestLess.getM();
}

long long Euler070(long long N)
{
    //8319823
    long long ri;
    long long rv;
    dec101_t rmin = N;
    long long idx = 0;
    dec101_t r;
    dec101_t p;
    std::vector<int> digita;
    std::vector<int> digitb;

    for(long long i = 2; i< N;i++)
    {
        p = phi(i);
        ri = (long long)p;
        rv = (long long)reverse_unumber(ri);
        digita = digits10(i);
        digitb = digits10(rv);
        std::sort(digita.begin(), digita.end());
        std::sort(digitb.begin(), digitb.end());
        if (digita == digitb)
        {
            r = dec101_t(i) / p;
            if (r < rmin)
            {
                rmin = r; idx=i;
                std::cout << i << " " << ri << " " << r << std::endl;
            }
        }
    }
    return idx;;
}

long long Euler069(long long N)
{
    //510510
    dec101_t r;
    dec101_t rmax = 0;
    long long idx = 0;
    for(long long i = 2; i<= N;i++)
    {
        r = i / phi(i);
        if (r > rmax)
        {
            rmax = r; idx=i;
        }
    }
    return  idx;;
}


void fillLine(  unsigned int N, unsigned int tripletSum, std::set<std::string> &result,
                unsigned int pos, std::vector<unsigned int> inner, std::vector<unsigned int> outer, unsigned int used)
{
  // inner ring completely filled, just one cell of the outer ring left
  if (pos == N - 1)
  {
    // check last line
    outer[N - 1] = tripletSum - (inner[0] + inner[N - 1]);
    unsigned int mask = 1 << outer[N - 1];
    if ((used & mask) != 0)
      return;

    // first element of outer ring must be the smallest
    for (auto x : outer)
      if (x < outer[0])
        return;

    // build string
    std::string id;
    for (unsigned int i = 0; i < N; i++)
      id += std::to_string(outer[i]) + std::to_string(inner[i]) + std::to_string(inner[(i + 1) % N]);

    // will be alphabetically ordered
    result.insert(id);
    return;
  }

  // move a number between 1 and 2*size into one of the inner cells of the n-gon
  for (unsigned int i = 1; i <= 2*N; i++)
  {
    // fill a cell of the inner ring
    unsigned int innerMask = 1 << i;
    // is that number still available ?
    if ((innerMask & used) != 0)
      continue;

    // occupy cell
    inner[pos + 1] = i;
    unsigned int nextUsed = used | innerMask;

    // compute the according cell in the outer ring
    outer[pos] = tripletSum - (inner[pos] + i);
    unsigned int outerMask = 1 << outer[pos];
    // is that number still available ?
    if ((nextUsed & outerMask) != 0)
      continue;
    nextUsed |= outerMask;

    // next line
    fillLine(N, tripletSum, result, pos + 1, inner, outer, nextUsed);
  }
}

std::string Euler068(long long N)
{
    //6531031914842725
    //https://euler.stephan-brumme.com/68/
    N = 5;
    unsigned int tripletSum = 14;
    std::set<std::string> result;

    // generate the inner and outer ring
    std::vector<unsigned int> inner(N);
    std::vector<unsigned int> outer(N);
    // a triplet consists of inner[a], inner[(a+1) % (2*a)], outer[a]

    // generate a bitmask of allowed numbers (0 = still available, 1 = already used / disallowed)
    unsigned int allowed = 0;
    for (unsigned int i = 1; i <= 2 * N; i++)
    allowed |= 1 << i;
    allowed = ~allowed;

    // fill first cell of inner ring
    for (unsigned int i = 1; i <= 2*N; i++)
    {
        inner[0] = i;
        // fill remaining cells
        fillLine((unsigned int)N, tripletSum, result, 0, inner, outer, allowed | (1 << i));
    }

    std::string  s;
    for (auto r : result)
    {
        std::cout << r << std::endl;
        s = r;
    }

    return s;
}

void visit_triangle(long long N, long long i, long long j, mat::matrix& mmaxpaths, mat::matrix& mnumber, mat::matrix& mvisit);
long long Euler067(long long N)
{
    // 7273
    // see Euler018
    int NROWS = (int)N;
    mat::matrix mnumber(NROWS, NROWS, 0);
    mat::matrix mmaxpaths(NROWS, NROWS, 0);
    mat::matrix mvisit(NROWS, NROWS, 0);
    //59
    //73 41
    //52 40 09
    //26 53 06 34
    //10 51 87 86 81
    //61 95 66 57 25 68
    //90 81 80 38 92 67 73
    // ...
    uinteger_t sum = 0;
    std::ifstream is("p067_triangle.txt", std::ifstream::in);
    if (is)
    {
        // get length of file:
        is.seekg(0, is.end);
        int length = (int)is.tellg();
        is.seekg(0, is.beg);

        char* buffer = new char[length];

        //std::cout << "Reading " << length << " characters... " << std::endl;;
        // read data as a block:
        is.read(buffer, length);
        if (is)
        {
            long num;
            int pos = 0;
            for (int i = 0; i < NROWS; i++)
            {
                num = 0;
                for (int j = 0; j < i + 1; j++)
                {
                    num = 10 * (buffer[pos] - '0') + (buffer[pos + 1] - '0');
                    pos += 3;
                    mnumber(i, j) = num;
                    std::cout << i << " " << j << " " << num << std::endl;
                }
            }
        }
        is.close();
    }
    visit_triangle(NROWS, 0, 0, mmaxpaths, mnumber, mvisit);
    return mmaxpaths(0, 0);
}

long long Euler066(long long N)
{
    // Pell equation
    // 661
    // https://euler.stephan-brumme.com/66/
    unsigned int limit = (unsigned int)N;

    unsigned int bestD = 2;
    uinteger_t bestX = 3;

    for (unsigned int d = 3; d <= limit; d++)
    {
        unsigned int root = (unsigned int)sqrt(d);
        if (root * root == d) continue;

        // see problem 64
        unsigned int a = root;
        unsigned int numerator = 0;
        unsigned int denominator = 1;

        // keep only the most recent 3 numerators and denominators while diverging
        uinteger_t x[3] = { 0, 1, root }; // numerators
        uinteger_t y[3] = { 0, 0, 1 };    // denominators

        // find better approximations until the exact solution is found
        while (true)
        {
            numerator = denominator * a - numerator;
            denominator = (d - numerator * numerator) / denominator;
            a = (root + numerator) / denominator;

            // x_n = a * x_n_minus_1 + x_n_minus_2
            x[0] = std::move(x[1]);
            x[1] = std::move(x[2]);
            x[2] = x[1] * a + x[0];

            // y_n = a * y_n_minus_1 + y_n_minus_2
            y[0] = std::move(y[1]);
            y[1] = std::move(y[2]);
            y[2] = y[1] * a + y[0];

            // avoid subtraction (to keep BigNum's code short)
            // x*x - d*y*y = 1
            // x*y         = 1 + d*y*y
            auto leftSide = x[2] * x[2];
            auto rightSide = y[2] * y[2] * d + 1;

            // solved it
            if (leftSide == rightSide)
                break;
        }

        // biggest x so far ?
        if (bestX < x[2])
        {
            bestX = x[2];
            bestD = d;
        }
    }
    return (long long) bestD;
}

long long Euler065(long long N)
{
    // 272
    //https://euler.stephan-brumme.com/65/
    // to save memory we dont keep all numerators, only the latest three
    uinteger_t numerators[3] = {    0,   // dummy, will be overwritten immediately
                                    1,   // always 1
                                    2 }; // the first number of the continuous fraction ("before the semicolon")

    for (unsigned int index = 2; index <= N; index++)
    {
        // e = [2; 1,2,1, 1,4,1, ... 1,2k,1, ...]
        unsigned int fractionNumber = 1;
        if (index % 3 == 0)
            fractionNumber = (index / 3) * 2;

        numerators[0] = std::move(numerators[1]);
        numerators[1] = std::move(numerators[2]);
        if (fractionNumber == 1)
            numerators[2] = numerators[0] + numerators[1];
        else
            numerators[2] = numerators[0] + numerators[1] * fractionNumber;
    }

    uinteger_t sum = 0;
    auto x = numerators[2];
    while (x > 0)
    {
        sum += x % 10;
        x /= 10;
    }
    return (long long)sum;
}


//https://euler.stephan-brumme.com/64/
unsigned int getPeriodLength(unsigned int x)
{
    unsigned int root = (unsigned int)sqrt(x);
    if (root * root == x)
        return 0;

    unsigned int a = root;
    unsigned int numerator = 0;     // initially zero, e.g. 4 will appear in second iteration of sqrt(23)
    unsigned int denominator = 1;   // initially one,  e.g. 7 will appear in second iteration of sqrt(23)
    unsigned int period = 0;

    // terminate when we see the same triplet (a, numerator, denominator) a second time
    // to me it wasn't obvious that this happens exactly when a == 2 * root
    // but thanks to Wikipedia for that trick ...
    while (a != 2 * root)
    {
        numerator = denominator * a - numerator;
        // must be integer divisions !
        denominator = (x - numerator * numerator) / denominator;
        a = (root + numerator) / denominator;
        period++;
    }
    return period;
}

long long Euler064(long long N)
{
    //1322
    unsigned int numOdd = 0;
    for (unsigned int i = 2; i <= N; i++) // 0 and 1 are perfect squares
    {
        unsigned int period = getPeriodLength(i);
        if (period % 2 == 1)
            numOdd++;
    }
    return numOdd;
}

long long Euler063(long long N)
{
    //49
    long long cnt = 1;
    uinteger_t npow;
    long long ndigit;

    for (long long a = 2; a < N; a++)
    {
        for (long long b = 1; b < N; b++)
        {
            npow = upow(a, b);
            ndigit = udigits10(npow).size();
            if (ndigit == b)
            {
                cnt++;
                std::cout << npow << " " << cnt << std::endl;
            }
            if (ndigit > b) break;
        }
    }
    return cnt;
}

long long Euler062(long long N)
{
    struct vlist
    {
        std::vector<uinteger_t> v;
    };

    //127035954683
    uinteger_t cube;
    uinteger_t smallest = N * N * N;
    std::map<std::vector<int>, vlist> vmap3;

    std::vector<int> v;
    for (long long i = 1; i < N; i++)
    {
        cube = i * i * i;
        v = udigits10(cube);
        std::sort(v.begin(), v.end());
        if (vmap3.find(v) == vmap3.end())
        {
            vlist vl;
            vl.v.push_back(cube);
            vmap3[v] = vl;
        }
        else
        {
            vlist vl = vmap3[v];
            vl.v.push_back(cube);
            vmap3[v] = vl;
        }
    }

    for (auto& [a, b] : vmap3)
    {
        if (b.v.size() == 5)
        {
            for (size_t i = 0; i < b.v.size(); i++)
            {
                if (b.v[i] < smallest)
                    smallest = b.v[i];
            }
        }
    }
    return (long long)smallest;
}

long long cyclic(std::map<long long, bool>& vmap3, std::map<long long, bool>& vmap4,
            std::map<long long, bool>& vmap5, std::map<long long, bool>& vmap6,
            std::map<long long, bool>& vmap7, std::map<long long, bool>& vmap8)
{
    long long  sum = 0;
    long long t3;
    long long t4;
    long long t5;
    long long t6;
    long long t7;
    long long t8;

    for (long long i=10; i <= 99; i++)
    for (long long j=0; j <= 99; j++)
    {
        t3 = 100*i + j;
        if (vmap3.find(t3) != vmap3.end())
        {
            //std::cout << "t3  " << t3 << std::endl;
            for (long long k=0; k <= 99; k++)
            {
                t4 = 100*j + k;
                if (vmap4.find(t4) != vmap4.end())
                {
                    //std::cout << "t4  " << t4 << std::endl;
                    for (long long l=0; l <= 99; l++)
                    {
                        t5 = 100*k + l;
                        if (vmap5.find(t5) != vmap5.end())
                        {
                            //std::cout << "t5  " << t5 << std::endl;
                            for (long long m=0; m <= 99; m++)
                            {
                                t6 = 100*l + m;
                                if (vmap6.find(t6) != vmap6.end())
                                {
                                    //std::cout << "t6  " << t6 << std::endl;
                                    for (long long n=0; n <= 99; n++)
                                    {
                                        t7 = 100*m +n;
                                        if (vmap7.find(t7) != vmap7.end())
                                        {
                                            //std::cout << "t7  " << t7 << std::endl;
                                            for (long long o=0; o <= 99;o++)
                                            {
                                                t8 = 100*n + o;
                                                if (vmap8.find(t8) != vmap8.end())
                                                {
                                                    std::cout << t3 << " " << t4 << " " << t5 << " "<< t6 << " "<< t7 << " "<< t8 << std::endl;
                                                    if (o == i)
                                                    {
                                                        long long sum = t3+t4+t5+t6+t7+t8;
                                                        std::cout << "DONE " << sum << std::endl;
                                                        return sum;
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    return sum;
}

long long Euler061(long long N)
{
    long long sum = 0;
    long long t3;
    long long t4;
    long long t5;
    long long t6;
    long long t7;
    long long t8;
    std::map<long long, bool> vmap3;
    std::map<long long, bool> vmap4;
    std::map<long long, bool> vmap5;
    std::map<long long, bool> vmap6;
    std::map<long long, bool> vmap7;
    std::map<long long, bool> vmap8;

    for (long long i=1; i < N; i++)
    {
        t3 = i*(i+1)/2;
        if (t3 > 9999) break;
        if (t3 >= 1000)
        {
            if (vmap3.find(t3) == vmap3.end()) vmap3[t3] = true;
        }
    }
    for (long long i=1; i < N; i++)
    {
        t4 = i*i;
        if (t4 > 9999) break;
        if (t4 >= 1000)
        {
            if (vmap4.find(t4) == vmap4.end()) vmap4[t4] = true;
        }
    }
    for (long long i=1; i < N; i++)
    {
        t5 = i*(3*i-1)/2;
        if (t5 > 9999) break;
        if (t5 >= 1000)
        {
            if (vmap5.find(t5) == vmap5.end()) vmap5[t5] = true;
        }
    }
    for (long long i=1; i < N; i++)
    {
        t6 = i*(2*i-1);
        if (t6 > 9999) break;
        if (t6 >= 1000)
        {
            if (vmap6.find(t6) == vmap6.end()) vmap6[t6] = true;
        }
    }
    for (long long i=1; i < N; i++)
    {
        t7 = i*(5*i-3)/2;
        if (t7 > 9999) break;
        if (t7 >= 1000)
        {
            if (vmap7.find(t7) == vmap7.end()) vmap7[t7] = true;
        }
    }
    for (long long i=1; i < N; i++)
    {
        t8 = i*(3*i-2);
        if (t8 > 9999) break;
        if (t8 >= 1000)
        {
            if (vmap8.find(t8) == vmap8.end()) vmap8[t8] = true;
        }
    }

    // one order  other possible!
    sum = cyclic(vmap3,vmap4,vmap5,vmap6,vmap7,vmap8); if (sum > 0) return sum;
    sum = cyclic(vmap3,vmap4,vmap5,vmap6,vmap8,vmap7); if (sum > 0) return sum;
    sum = cyclic(vmap3,vmap4,vmap5,vmap7,vmap6,vmap8); if (sum > 0) return sum;
    sum = cyclic(vmap3,vmap4,vmap5,vmap7,vmap8,vmap6); if (sum > 0) return sum;
    sum = cyclic(vmap3,vmap4,vmap5,vmap8,vmap6,vmap7); if (sum > 0) return sum;
    sum = cyclic(vmap3,vmap4,vmap5,vmap8,vmap7,vmap6); if (sum > 0) return sum;
    sum = cyclic(vmap3,vmap4,vmap6,vmap5,vmap7,vmap6); if (sum > 0) return sum;
    sum = cyclic(vmap3,vmap4,vmap6,vmap5,vmap8,vmap7); if (sum > 0) return sum;
    sum = cyclic(vmap3,vmap4,vmap6,vmap7,vmap5,vmap8); if (sum > 0) return sum;
    sum = cyclic(vmap3,vmap4,vmap6,vmap7,vmap8,vmap5); if (sum > 0) return sum;
    sum = cyclic(vmap3,vmap4,vmap6,vmap8,vmap5,vmap7); if (sum > 0) return sum;
    sum = cyclic(vmap3,vmap4,vmap6,vmap8,vmap7,vmap5); if (sum > 0) return sum;
    sum = cyclic(vmap3,vmap4,vmap7,vmap5,vmap6,vmap8); if (sum > 0) return sum;
    sum = cyclic(vmap3,vmap4,vmap7,vmap5,vmap8,vmap6); if (sum > 0) return sum;
    sum = cyclic(vmap3,vmap4,vmap7,vmap6,vmap5,vmap8); if (sum > 0) return sum;
    sum = cyclic(vmap3,vmap4,vmap7,vmap6,vmap8,vmap5); if (sum > 0) return sum;
    sum = cyclic(vmap3,vmap4,vmap7,vmap8,vmap5,vmap6); if (sum > 0) return sum;
    sum = cyclic(vmap3,vmap4,vmap7,vmap8,vmap6,vmap5); if (sum > 0) return sum;

    return 0;
}


long long concate(long long  a, long long b)
{
    long long r = a;
    long long t = b;
    while (t > 0)
    {
        r = 10*r;
        t = t/10;
    }
    r += b;
    return r;
}

long long Euler060(long long N)
{
    std::map<std::vector<long long>, long long> vamp2set;
    std::map<std::vector<long long>, long long> vamp3set;
    std::map<std::vector<long long>, long long> vamp4set;
    std::map<std::vector<long long>, long long> vamp5set;

    long long  sum;
    long long  cnt = 0;
    long long  minsum = 2*N;
    std::vector<long long> v;

    for (int a=3; a < N; a+=2)
    for (int b=a+2; b < N; b+=2)
    {
        cnt++;
        if (  (a%5 != 0) && (b%5 != 0) &&
              is_prime(a) && is_prime(b)  &&
              is_prime(concate(a, b)) && is_prime(concate(b, a))
            )
        {
            sum = a+b;
            v = {a,b};
            vamp2set[v] = sum;
            std::cout << "sum2  " << sum << " " <<  a << " " <<  b << std::endl;

            if (sum < minsum)
            {
                minsum = sum;
                std::cout << "minsum2  " << minsum << " " <<  a << " " <<  b << std::endl;
            }
        }
        if (cnt%100001==0) std::cout << "cnt  " << cnt << " " <<  a << " " <<  b<< std::endl;
    }

    minsum = 3*N;
    for(auto& [v, s] : vamp2set)
    {
        long long a = v[0]; long long b=v[1];
        for (long long c=b+2; c < N; c+=2)
        {
            cnt++;
            if (  (c%5 != 0) && is_prime(c)  &&
                  is_prime(concate(a, c)) && is_prime(concate(c, a)) &&
                  is_prime(concate(b, c)) && is_prime(concate(c, b))
                )
            {
                sum = a+b+c;
                std::vector <long long> v3  = {a, b, c};
                vamp3set[v3] = sum;
                std::cout << "sum3  " << sum << " " <<  a << " " <<  b << " " <<  c << std::endl;

                if (sum < minsum)
                {
                    minsum = sum;
                    std::cout << "minsum3 " << minsum << " " <<  a << " " <<  b << " " <<  c <<   std::endl;
                }
            }
            if (cnt%100001==0) std::cout << "cnt  " << cnt << " " <<  a << " " <<  b << " " <<  c <<  std::endl;
        }
    }

    minsum = 4*N;
    for(auto& [v, s] : vamp3set)
    {
        long long a = v[0];
        long long b=v[1];
        long long c = v[2];

        for (long long d=c+2; d < N; d+=2)
        {
            if (  (d%5 != 0) && is_prime(d) )
            {
                cnt++;
                if (
                      is_prime(concate(a, d)) && is_prime(concate(d, a))  &&
                      is_prime(concate(b, d)) && is_prime(concate(d, b))  &&
                      is_prime(concate(c, d)) && is_prime(concate(d, c))
                    )
                {
                    sum = a+b+c+d;
                    std::vector <long long> v4 = {a, b, c, d};
                    vamp4set[v4] = sum;
                    if (sum < minsum)
                    {
                        minsum = sum;
                        std::cout << "minsum4 " << minsum << " " <<  a << " " <<  b << " " <<  c << " " << d <<  std::endl;
                    }
                }
            }
        }
    }

    minsum = 5*N;
    for(auto& [v, s] : vamp4set)
    {
        long long a = v[0];
        long long b = v[1];
        long long c = v[2];
        long long d = v[3];

        for (long long e=d+2; e < N; e+=2)
        {
            if (  (e%5 != 0) && is_prime(e) )
            {
                cnt++;
                if (
                      is_prime(concate(a, e)) && is_prime(concate(e, a))  &&
                      is_prime(concate(b, e)) && is_prime(concate(e, b))  &&
                      is_prime(concate(c, e)) && is_prime(concate(e, c))  &&
                      is_prime(concate(d, e)) && is_prime(concate(e, d))
                    )
                {
                    sum = a+b+c+d+e;
                    std::vector <long long> v5 = {a, b, c, d, e};
                    vamp5set[v5] = sum;
                    if (sum < minsum)
                    {
                        minsum = sum;
                        std::cout << "minsum5 " << minsum << " " <<  a << " " <<  b << " " <<  c << " " << d << " " << e  << std::endl;
                    }
                }
            }
        }
    }
    //minsum5 26033 13 5197 5701 6733 8389
    return minsum;
}


long long Euler059()
{
    double sc = 0;
    int num;
    std::vector<int> bestkey;
    std::vector<int> vencrypted;
    std::map<int, long long> vmap_decrypted;

    std::ifstream is("p059_cipher.txt", std::ifstream::in);
    if (is)
    {
        // get length of file:
        is.seekg (0, is.end);
        int length = (int)is.tellg();
        is.seekg (0, is.beg);

        char * buffer = new char [length];

        //std::cout << "Reading " << length << " characters... "<< std::endl;;
        is.read (buffer,length);
        std::string s;
        std::vector<std::string> v;

        if (is)
        {
            //36,22,80,0,0,4,23,25,19,17,88,4,4,19,21,
            for (int j=0; j < length; j++)
            {
                if (buffer[j]==',')
                {
                    if (s.size()>0)
                    {
                        v.push_back(s);
                        s.clear();
                    }
                }
                else if (buffer[j]!=',') s=s+buffer[j];
            }
            if (s.size()>0) v.push_back(s);
            std::cout  << v.size() << std::endl;
        }
        is.close();

        for (size_t i=0; i < v.size(); i++)
        {
            num = stoi(v[i]);
            vencrypted.push_back(num);
        }

        std::vector<int> letter {'e', 't', 'a', 'o', 'i', 'E', 'T', 'A', 'O', 'I', };
        std::vector<double> freq {0.1202, 0.910, 0.812, 0.768, 0.731, 0.1202, 0.910, 0.812, 0.768, 0.731};

        std::vector<int> key(v.size(), 0);
        double best_sc = 9999999;

        for (int a=0; a < 26; a++)
        for (int b=0; b < 26; b++)
        for (int c=0; c < 26; c++)
        {
            key.clear();
            vmap_decrypted.clear();

            while(key.size() < vencrypted.size() )
            {
                if (key.size() < vencrypted.size()) key.push_back('a'+a);
                if (key.size() < vencrypted.size()) key.push_back('a'+b);
                if (key.size() < vencrypted.size()) key.push_back('a'+c);
            }
            //s = {key[0], key[1], key[2]};

            for (size_t i=0; i < vencrypted.size(); i++)
            {
                num = vencrypted[i];
                num ^= key[i];
                vmap_decrypted[num]++;
            }

            sc = 0;
            for (size_t i=0; i < 10; i++)
            {
                double r = (double)vmap_decrypted[letter[i]]; r = r / ((double)vencrypted.size());
                sc +=  std::abs(r-freq[i]) * std::abs(r-freq[i]);
            }
            if (sc < best_sc)
            {
                best_sc = sc;
                bestkey = key;
                std::cout << "best  " << s << " "<< best_sc << " " <<  std::endl;

                for (size_t i=0; i < vencrypted.size(); i++)
                {
                    num = vencrypted[i];
                    num ^= key[i];
                    //if (i == 0) std::cout << s << " ";
                    if (i < 80) std::cout << (char)num;
                }
                std::cout <<  std::endl;
            }
        }
    }

    long long sum = 0;
    for (size_t i=0; i < vencrypted.size(); i++)
    {
        num = vencrypted[i] ^ bestkey[i];
        sum += num;
    }
    return sum;
}

void spiral_corner_l(long long n, std::vector<long long>& v)
{
    if (n==1){v.push_back(1); return;}
    long long delta = n-1;
    v.push_back(n*n - 1 * delta);
    v.push_back(n*n - 3 * delta);
}
void spiral_corner_r(long long n, std::vector<long long>& v)
{
    if (n==1){v.push_back(1); return;}
    long long delta = n-1;
    v.push_back(n*n);
    v.push_back(n*n - 2 * delta);
}

long long Euler058(long long N)
{
    std::vector<long long> vr;
    std::vector<long long> vl;
    long long cnt_p = 0;
    long long cnt = 0;
    long long last_cnt = 0;
    for (long long i = 1; i < N; i+=2)
    {
        spiral_corner_l(i, vl);
        spiral_corner_r(i, vr);
        for (size_t j = last_cnt; j < vl.size(); j++)
            if (is_prime(vl[j])) cnt_p++;
        for (size_t j = last_cnt; j < vr.size(); j++)
            if (is_prime(vr[j])) cnt_p++;
        last_cnt = vl.size();
        cnt = vl.size() + vr.size() - 1;
        double r = ((double)cnt_p) / ((double)cnt) ;

        std::cout << i << " " << r <<  std::endl;

        if (i> 2)
        if (r < 0.100) return i;
    }
    return 0;
}

dec101_t Euler826(long long N )
{
    //https://www.cut-the-knot.org/Curriculum/Probability/BOW5.shtml
    // 2/n+1 + (n-3)(7/18*(n+1)
    dec101_t cnt = 0;
    dec101_t r;
    dec101_t rn;
    dec101_t tdist = 0;
    for (long long i = 3; i < N; i++)
    {
        if (is_prime(i))
        {
            cnt++;
            rn = i;
            r = 2.0/(rn+1);
            r += (rn-3)*7/(18*(rn+1.0));
            tdist +=r;
        }
    }
    std::cout << " " << std::setprecision(20) << tdist/cnt << std::endl; //0.3889014797
    return tdist/cnt;
}

// TODO uint2048_t not enough large maybe
long long Euler057(long long N)
{
    // 153
    long long cnt = 0;
    RationalNumber ONE(BigIntegerONE);
    RationalNumber TWO(BigInteger(2));
    RationalNumber half = ONE/TWO;
    RationalNumber last_fraction = half;
    RationalNumber r;

    std::cout << 0 << " " << half << " " << cnt << std::endl;
    for (long long i = 0; i < N; i++)
    {
        r = ONE + ONE / (TWO + last_fraction);
        last_fraction = ONE / (TWO + last_fraction);

        auto va = udigits10(r.getN().toLongLong());
        auto vb = udigits10(r.getM().toLongLong());
        if (va.size() < vb.size())
        {
            cnt++;
            std::cout << i << " " << r << " " << cnt << std::endl;
        }
    }
    return cnt;
}


long long Euler056(long long N)
{
    uinteger_t sum = 0;
    uinteger_t maxsum = 0;
    uinteger_t c;
    std::vector<int> v;
    for (long long a = 1; a < N; a++)
    for (long long b = 1; b < N; b++)
    {
        c = upow(a, b);
        v = udigits10(c);
        sum = 0;
        for (size_t i = 0; i < v.size(); i++)
            sum += v[i];
        if (sum > maxsum) maxsum = sum;
    }
    return (long long)maxsum;
}

bool is_lych(long long k)
{
    uinteger_t sum = k;
    uinteger_t r;

    for (long long i = 1; i < 50; i++)
    {
        r = reverse_unumber(sum);
        sum+=r;
        if (is_upalindrome(sum)) return false;
    }
    return true;
}

long long Euler055(long long N)
{
    long long cnt = 0;
    for (long long r = 1; r < N; r++)
    {
        if (is_lych(r))
            cnt++;
    }
    return cnt;
}

long long Euler053(long long N)
{
    uinteger_t c;
    long long cnt = 0;

    for (long long n = 1; n <= N; n++)
    for (long long r = 1; r <= n; r++)
    {
        c = ucombination(r, n);
        if (c > 1000000)
            cnt++;
    }
    return cnt;
}

long long Euler052(long long N)
{
    bool ok;
    std::map<int, int> vmap_digit[6];

    for (long long i = 2; i <= N; i++)
    {
        ok = true;
        for (int j = 0; j < 6; j++)
            vmap_digit[j] = digitsmap10((j+1) * i);

        for (int j = 1; j < 6; j++)
        {
            if (vmap_digit[j] != vmap_digit[0])
            {
                ok = false;
                break;
            }
        }
        if (ok) return i;
    }
    return 0;
}

long long replace_digit(long long i, int digits, int j)
{
    long long t = 0;
    auto va = digits10(i);
    for (size_t i = 0; i < va.size(); i++)
    {
        if (va[i] == digits) va[i]=j;
    }
    for (size_t i = 0; i < va.size(); i++)
    {
        t = 10 * t + va[va.size() - i - 1];
    }
    return t;
}

long long Euler051(long long N, int CNT)
{
    long cnt=0;
    std::map<int, bool> vmap_digit;

    for (long long i = 2; i <= N; i++)
    {
        if (is_prime(i))
        {
            auto va = digits10(i);
            vmap_digit.clear();
            for (size_t i = 0; i < va.size(); i++)
            {
                vmap_digit[va[i]] = true;
            }

            for (int digits = 0; digits <= 9; digits++)
            {
                if (vmap_digit.find(digits) != vmap_digit.end())
                {
                    cnt = 0;
                    for (int j = 0; j <= 9; j++)
                    {
                        auto t = replace_digit(i, digits, j);
                        auto vt = digits10(t);
                        if (vt.size() == va.size())
                        {
                            if (is_prime(t))
                            {
                                cnt++;
                            }
                        }
                    }
                }
                if (cnt==CNT) return i;
            }
        }
    }
    return 0;
}

long long Euler050(long long N)
{
    long long cnt = 0;
    long long sum = 0;
    long long maxs = 0;
    long long maxcnt = 0;
    long long p ;
    for (long long i = 2; i < N; i++)
    {
        std::cout << i << " " << std::endl;
        sum = 0;
        cnt = 0;
        if (is_prime(i))
        {
            sum+=i;
            cnt++;
            p = i;
            while(true)
            {
                p = next_prime(p, N);
                if (p==2) break;
                sum+=p;
                cnt++;
                if (sum >= N) break;

                if (is_prime(sum))
                if (cnt > maxcnt)
                {
                    maxcnt = cnt;
                    maxs = sum;
                    std::cout << i << " maxcnt " << maxcnt << " " << sum << std::endl;
                }
            }
        }
    }
    return maxs;
}

bool is_permutation(long long a, long long b)
{
    auto va = digits10(a); // std::vector<int>
    auto vb = digits10(b);
    std::map<long long, long long> vmapa;
    std::map<long long, long long> vmapb;
    for (size_t i = 0; i < va.size(); i++)
    {
        vmapa[ va[i] ]++;
    }
    for (size_t i = 0; i < vb.size(); i++)
    {
        vmapb[ vb[i] ]++;
    }
    if (va.size()!=vb.size()) return false;
    if (vmapa.size()!=vmapb.size()) return false;
    for (auto& [aa, ab] : vmapa)
    {
        if (vmapb.find(aa) == vmapb.end()) return false;
    }
    return true;
}

long long Euler049(long long N)
{
    uinteger_t t = 0;
    uinteger_t sum = 0;
    for (long long i = 1000; i < N; i++)
    {
        std::cout << i << " " << std::endl;
        if (i!=1487)
        for (long long delta = 1; delta <= N/3; delta++)
        {
            if (i+2*delta < N)
            if (is_prime(i))
            if (is_prime(i+delta))
            if (is_prime(i+2*delta))
            if (is_permutation(i, i+delta))
            if (is_permutation(i, i+2*delta))
            {
                return i*100000000 + (i+delta)*10000 + i+2*delta;
            }
        }
    }
    return 0;
}

long long Euler048(long long N)
{
    uinteger_t t = 0;
    uinteger_t sum = 0;
    int cnt = 0;
    for (long long i = 1; i <= N; i++)
        t += upow(i, i);

    while(cnt < 10)
    {
        uinteger_t digit = t % 10;
        sum = sum + digit * upow(10, cnt);
        t = t - digit;
        t = t / 10;
        cnt++;
    }
    return (long long)sum;
}

long long Euler047(long long N)
{
    //134043 SLOW
    long long r[4];
    std::unordered_map<long long, long long> vmap_factors;
    std::unordered_map<long long, std::unordered_map<long long, long long>> vmap_cache;

    bool ok;
    for (long long i = 16; i <= N-4; i++)
    {
        ok = true;
        r[0] = i;
        r[1] = i+1;
        r[2] = i+2;
        r[3] = i+3;
        for (long long m = 0; m < 4; m++)
        {
            if (vmap_cache.find(i+m)!=vmap_cache.end())
            {
                vmap_factors = vmap_cache[i+m];
            }
            else
            {
                vmap_factors = unique_prime_factors(r[m], 4);
                vmap_cache[i+m] = vmap_factors ;
            }
            if (vmap_factors.size() != 4) ok=false;;
            if (ok==false) break;
        }
        if (i % 100==0) std::cout << i << " " << std::endl;
        if (ok) return i;
    }
    return 0;
}

long long Euler046(long long N)
{
    long long r;
    long long sq=0;
    for (long long i = 9; i <= N; i++)
    {
        if (i % 2 == 1)
        if (is_prime(i)==false)
        {
            for (long long j = 1; j <= N; j++)
            {
                sq = j*j;
                if (sq >= i) return i;
                r = i - 2*sq;
                if (r == 0) return i;
                if (is_prime(r)==true)
                {
                    break;
                }
            }
        }
    }
    return 0;
}

long long Euler045(long long N)
{
    //1533776805
    std::map<long long, long long> vmap3;
    std::map<long long, long long> vmap5;
    std::map<long long, long long> vmap8;
    std::map<long long, long long> vmap3_inv;
    std::map<long long, long long> vmap5_inv;
    std::map<long long, long long> vmap8_inv;

    long long r;
    long long k=0;
    for (long long i = 1; i <= N; i++)
    {
        vmap3[i] = (i*(i+1)/2);
        vmap3_inv[(i*(i+1)/2)] = i;
        vmap5[i] = (i*(3*i-1)/2);
        vmap5_inv[(i*(3*i-1)/2)] = i;
        vmap8[i] = (i*(2*i-1));
        vmap8_inv[(i*(2*i-1))] = i;
    }

    for (long long i = 286; i <= N; i++)
    {
        r = vmap3[i];
        if (vmap5_inv.find(r) != vmap5_inv.end())
        {
            if (vmap8_inv.find(r) != vmap8_inv.end())
            {
                std::cout << r << std::endl;
                k = r;
                break;
            }
        }

    }
    return  k;
}

long long Euler044(long long N)
{
    std::map<long long, long long> vmap;
    std::map<long long, long long> vmap_inv;
    long long minp = 10000000000;
    long long r;long long rr;
    for (long long i = 1; i <= N; i++)
    {
        vmap[i] = (i*(3*i-1)/2);
        vmap_inv[(i*(3*i-1)/2)] = i;
    }

    for (long long i = 1; i <= N; i++)
    {
        for (long long j = i+1; j <= N; j++)
        {
            r = vmap[j] - vmap[i];
            rr = vmap[j] + vmap[i];
            if (r < 0) r = -1*r;
            if (r > minp) break;

            std::cout << i << " " << r << " " << minp << std::endl;

            if (vmap_inv.find(r) != vmap_inv.end())
            {
                if (vmap_inv.find(rr) != vmap_inv.end())
                {
                 std::cout << r << std::endl;
                    if (r < minp) minp = r;
                }
            }
        }
    }
    return minp; // 5482660
}

long long Euler043()
{
    long long sum=0;
    long long t;

    {
        for (size_t i = 0; i <= 9; i++)
        for (size_t j = 0; j <= 9; j++)
        for (size_t k = 0; k <= 9; k++)
        {
            if ( (i!=j) && (i!=k) && (j!=k) )
            {
                for (size_t l = 0; l <= 9; l++)
                {
                    if ( (l!=i) && (l!=j) && (l!=k) && (((100*j +10*k + l) % 2) == 0) )
                    {
                        for (size_t m = 0; m <= 9; m++)
                        {
                            if ( (m!=i) && (m!=j) && (m!=k) && (m!=l) && (((100*k +10*l + m) % 3) == 0) )
                            {
                                for (size_t n = 0; n <= 9; n++)
                                {
                                    if ( (n!=i) && (n!=j) && (n!=k) && (n!=l) && (n!=m) && (((100*l +10*m + n) % 5) == 0) )
                                    {
                                        for (size_t o = 0; o <= 9; o++)
                                        {
                                            if ( (o!=i) && (o!=j) && (o!=k) && (o!=l) && (o!=m) && (o!=n) && (((100*m +10*n + o) % 7) == 0) )
                                            {
                                                for (size_t p = 0; p <= 9; p++)
                                                {
                                                    if ( (p!=i) && (p!=j) && (p!=k) && (p!=l) && (p!=m) && (p!=n) && (p!=o) && (((100*n +10*o + p) % 11) == 0) )
                                                    {
                                                        for (size_t q = 0; q <= 9; q++)
                                                        {
                                                            if ( (q!=i) && (q!=j) && (q!=k) && (q!=l) && (q!=m) && (q!=n) && (q!=o) && (q!=p) && (((100*o +10*p + q) % 13) == 0) )
                                                            {
                                                                for (size_t r = 0; r <= 9; r++)
                                                                {
                                                                    if ( (r!=i) && (r!=j) && (r!=k) && (r!=l) && (r!=m) && (r!=n) && (r!=o) && (r!=p) && (r!=q) && (((100*p +10*q + r) % 17) == 0) )
                                                                    {
                                                                        t = r + 10*q + 100*p + 1000*o  +  10000*n  +  100000*m  + 1000000*l  + 10000000*k  + 100000000*j  + 1000000000*i;
                                                                        {
                                                                            std::cout << t << std::endl;
                                                                            sum+= t;
                                                                        }
                                                                    }
                                                                }
                                                            }
                                                        }
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    return sum;
}


bool is_triangle(long long r)
{
    for (long long j = 1; j < 1000; j++)
    {
        if (2*r == j*(j+1)) return true;
    }
    return false;
}

long long Euler042()
{
    long long cnt = 0;
    std::ifstream is("p042_words.txt", std::ifstream::in);
    if (is)
    {
        // get length of file:
        is.seekg (0, is.end);
        int length = (int)is.tellg();
        is.seekg (0, is.beg);

        char * buffer = new char [length];

        //std::cout << "Reading " << length << " characters... "<< std::endl;;
        is.read (buffer,length);
        std::string s;
        std::vector<std::string> v;

        if (is)
        {
            //"A","ABILITY","ABLE","A
            for (int j=0; j < length; j++)
            {
                if (buffer[j]=='\"')
                {
                    if (s.size()>0)
                    {
                        v.push_back(s);
                        s.clear();
                    }
                }
                else if (buffer[j]!=',') s=s+buffer[j];
            }
            std::cout  << v.size() << std::endl;
        }
        is.close();

        std::sort(v.begin(), v.end());
        int sc = 0;
        for (size_t i=0; i < v.size(); i++)
        {
            s = v[i];
            sc = 0;
            for (size_t j=0; j < s.size(); j++)
            {
                sc += (1 + s[j] - 'A');
            }
            if (is_triangle(sc))
                cnt++;
        }
    }
    return cnt;
}

bool is_pandigital(long long r)
{
    std::vector<int> v;
    std::map<int, int> vm;

    v = digits10(r);
    for (size_t j = 0; j < v.size(); j++)
    {
        if (v[j] == 0) return false;
        vm[v[j]-1]++;
        if (vm[v[j]-1] > 1) return false;
    }
    for (size_t j = 0; j < vm.size(); j++)
    {
        if (vm[(int)j] != (int)1) return false;
    }
    return true;
}

long long Euler041(long long N)
{
    long long maxn=0;
    for (long long t = 1; t<=N; t++)
    {
        if (is_pandigital(t))
        if (is_prime(t))
        {
            std::cout << t << std::endl;
            if (t > maxn)  maxn = t;
        }
    }
    return maxn;
}


long long d_k(long long k)
{
    long long b;
    std::vector<int> v = {  9,
                            9 + 2*90,
                            9 + 2*90 + 3*(1000-100 + 0),
                            9 + 2*90 + 3*(1000-100 + 0) + 4*(10000-1000 + 0),
                            9 + 2*90 + 3*(1000-100 + 0) + 4*(10000-1000 + 0) + 5*(100000-10000 + 0),
                            9 + 2*90 + 3*(1000-100 + 0) + 4*(10000-1000 + 0) + 5*(100000-10000 + 0) + 6*(1000000-100000 + 0)
    };
    //0.123456789 101112131415161718192021...99
    //          9                             9+2*90
    if (k<=v[0]) return k;
    if ((k>v[0]) && (k<= v[1]))
    {
        b = v[0];
        if ((k-b-1) %2 == 0)    return digit_at(0, 10 + (k-b-1) / 2);
        else                    return digit_at(1, 10 + (k-b-1) / 2);
    }
    if ((k>v[1]) && (k<v[2]))
    {
        b = v[1];
        if ((k-b-1) %3 == 0)       return digit_at(0, 100 + (k-b-1) / 3);
        else if ((k-b-1) %3 == 1)  return digit_at(1, 100 + (k-b-1) / 3);
        else                       return digit_at(2, 100 + (k-b-1) / 3);
    }
    if ((k>v[2]) && (k<v[3]))
    {
        b = v[2];
        if ((k-b-1) %4 == 0)       return digit_at(0, 1000 + (k-b-1) / 4);
        else if ((k-b-1) %4 == 1)  return digit_at(1, 1000 + (k-b-1) / 4);
        else if ((k-b-1) %4 == 2)  return digit_at(2, 1000 + (k-b-1) / 4);
        else                       return digit_at(3, 1000 + (k-b-1) / 4);
    }
    if ((k>v[3]) && (k<v[4]))
    {
        b = v[3];
        if ((k-b-1) %5 == 0)      return digit_at(0, 10000 + (k-b-1) / 5);
        else if ((k-b-1) %5 == 1) return digit_at(1, 10000 + (k-b-1) / 5);
        else if ((k-b-1) %5 == 2) return digit_at(2, 10000 + (k-b-1) / 5);
        else if ((k-b-1) %5 == 3) return digit_at(3, 10000 + (k-b-1) / 5);
        else                      return digit_at(4, 10000 + (k-b-1) / 5);
    }
    if ((k>v[4]) && (k<v[5]))
    {
        b = v[4];
        if ((k-b-1) %6 == 0)      return digit_at(0, 100000 + (k-b-1) / 6);
        else if ((k-b-1) %6 == 1) return digit_at(1, 100000 + (k-b-1) / 6);
        else if ((k-b-1) %6 == 2) return digit_at(2, 100000 + (k-b-1) / 6);
        else if ((k-b-1) %6 == 3) return digit_at(3, 100000 + (k-b-1) / 6);
        else if ((k-b-1) %6 == 4) return digit_at(4, 100000 + (k-b-1) / 6);
        else                      return digit_at(5, 100000 + (k-b-1) / 6);
    }
    if (k>v[5] )
    {
        b = v[5];
        if ((k-b) %7 == 0)        return digit_at(0, 1000000 + (k-b-1) / 7);
        else if ((k-b-1) %7 == 1) return digit_at(1, 1000000 + (k-b-1) / 7);
        else if ((k-b-1) %7 == 2) return digit_at(2, 1000000 + (k-b-1) / 7);
        else if ((k-b-1) %7 == 3) return digit_at(3, 1000000 + (k-b-1) / 7);
        else if ((k-b-1) %7 == 4) return digit_at(4, 1000000 + (k-b-1) / 7);
        else if ((k-b-1) %7 == 5) return digit_at(5, 1000000 + (k-b-1) / 7);
        else                      return digit_at(6, 1000000 + (k-b-1) / 7);
    }
    return 0;
}

long long Euler040()
{
    std::cout << d_k(1) << std::endl;
    std::cout << d_k(10) << std::endl;
    std::cout << d_k(100) << std::endl;
    std::cout << d_k(1000) << std::endl;
    std::cout << d_k(10000) << std::endl;
    std::cout << d_k(100000) << std::endl;
    std::cout << d_k(1000000) << std::endl;
    return d_k(1) * d_k(10) * d_k(100) * d_k(1000) * d_k(10000) * d_k(100000) * d_k(1000000);
}

long long Euler039(long long N)
{
    long long maxn=0;
    long long cnt = 0;
    long long pmax = 0;
    long long r;
    for (long long p = 3; p < N; p++)
    {
        cnt = 0;
        for (long long a = 1; a < N; a++)
        {
            for (long long b = 1; b < N; b++)
            {
                long long c = a*a + b*b;
                r = (long long)std::sqrt(c); // trunc
                if (r*r == c)
                {
                    if (a+b+r == p)
                    {
                        cnt++;
                    }
                }
            }
        }
        if (cnt > maxn)
        {
            std::cout << p << std::endl;
            maxn = cnt; pmax = p;
        }
    }
    return pmax;
}

long long max_power10(long long k)
{
    std::vector<int> v;
    long long r = k;
    int digit;
    while (r > 0)
    {
        digit = r % 10;
        v.push_back(digit);
        r = r/10;
    }
    return v.size();
}
bool is_pandigital_9(long long k, long long& r)
{
    r = 0;
    long long maxpower10;
    std::vector<int> v;
    std::map<int, int> vm;
    for (long long i = 1; i < 10; i++)
    {
        maxpower10 = max_power10(k * i);
        r = r * (long long)upow(10, maxpower10);
        r += k * i;
        if (r > 987654321) return false;
        v = digits10(r);
        if (v.size() > 9) return false;
        vm.clear();
        for (size_t j = 0; j < v.size(); j++)
        {
            if (v[j] == 0) return false;
            vm[v[j]-1]++;
            if (vm[v[j]-1] > 1) return false;
        }
        if (v.size() == 9) break;
    }
    if (r < 123456789) return false;
    return true;
}

long long Euler038(long long N)
{
    long long maxn=0;
    long long r;
    for (long long i = 2; i < N; i++)
    {
        if (is_pandigital_9(i, r))
        {
            std::cout << r << std::endl;
            if (r > maxn)  maxn = r;
        }
    }
    return maxn;
}



long long Euler037(long long N)
{
    long long sum=0;
    for (long long i = 11; i < N; i++)
    {
        if (is_truncatable_prime(i))
        {
            std::cout << i << std::endl;
            sum += i;
        }
    }
    return sum;
}

long long Euler036(long long N)
{
    long long sum=0;
    for (long long i = 1; i < N; i++)
    {
        if ((i%10 != 0) && (i%2 != 0))
        if (is_palindrome_base10(i) && is_palindrome_base2(i))
        {
            std::cout << i << std::endl;
            sum += i;
        }
    }
    return sum;
}

long long Euler035(long long N)
{
    long long r;
    long long r1;long long r2;long long r3;long long r4;
    long long d1;
    long long d2;
    long long d3;long long d4;long long d5;long long d6;
    std::map<long long , bool> vmap;

    for(long long p=2;p<=N;p++)
    {
        if (is_prime(p))
        {
            if ((p>=2) && (p<10) )
            {
                vmap[p] = true;
                std::cout << p << std::endl;
            }

            if ((p>=10) && (p<100) )
            {
                d1 = p / 10;
                d2 = p % 10;
                r = 10*d2 + d1;
                if (is_prime(r))
                {
                    vmap[p] = true;
                    vmap[r] = true;
                    std::cout << p << std::endl;
                }
            }

            if ((p>=100) && (p<1000) )
            {
                d1 = p / 100;
                d2 = p - d1 * 100; d2 = d2 / 10;
                d3 = p - d1 * 100 - 10 * d2;
                r = 100*d2 + 10*d3 + d1;
                if (is_prime(r))
                {
                    r1 = 100*d3 + 10*d1 + d2;
                    if (is_prime(r1))
                    {
                        vmap[p] = true;
                        vmap[r] = true;
                        vmap[r1] = true;
                        std::cout << p << std::endl;
                    }
                }
            }

            if ((p>=1000) && (p<10000))
            {
                // 1234
                d1 = p / 1000;
                d2 = p - d1 * 1000; d2 = d2 / 100;
                d3 = p - d1 * 1000 - d2 * 100;  d3 = d3 / 10;
                d4 = p - d1 * 1000 - d2 * 100 - d3*10;

                // 4123
                // 3412
                // 2341
                r  = 1000*d4 + 100*d1 + 10*d2 + d3;
                r1 = 1000*d3 + 100*d4 + 10*d1 + d2;
                r2 = 1000*d2 + 100*d3 + 10*d4 + d1;
                if (is_prime(r))
                if (is_prime(r1))
                if (is_prime(r2))
                {
                    vmap[p] = true;
                    vmap[r] = true;
                    vmap[r1] = true;
                    vmap[r2] = true;
                    std::cout << p << std::endl;
                }
            }

            if ((p>=10000) && (p<100000))
            {
                // 12345
                d1 = p / 10000;
                d2 = p - d1 * 10000; d2 = d2 / 1000;
                d3 = p - d1 * 10000 - d2 * 1000;  d3 = d3 / 100;
                d4 = p - d1 * 10000 - d2 * 1000 - d3 * 100;  d4 = d4 / 10;
                d5 = p - d1 * 10000 - d2 * 1000 - d3 * 100 - d4 * 10;

                // 51234
                // 45123
                // 34512
                // 23451
                r  = 10000*d5 + 1000*d1 + 100*d2 + 10*d3 + d4;
                r1 = 10000*d4 + 1000*d5 + 100*d1 + 10*d2 + d3;
                r2 = 10000*d3 + 1000*d4 + 100*d5 + 10*d1 + d2;
                r3 = 10000*d2 + 1000*d3 + 100*d4 + 10*d5 + d1;
                if (is_prime(r))
                if (is_prime(r1))
                if (is_prime(r2))
                if (is_prime(r3))
                {
                    vmap[p] = true;
                    vmap[r] = true;
                    vmap[r1] = true;
                    vmap[r2] = true;
                    vmap[r3] = true;
                    std::cout << p << std::endl;
                }

            }

            if ((p>=100000) && (p<1000000))
            {
                // 123456
                d1 = p / 100000;
                d2 = p - d1 * 100000; d2 = d2 / 10000;
                d3 = p - d1 * 100000 - d2 * 10000;  d3 = d3 / 1000;
                d4 = p - d1 * 100000 - d2 * 10000 - d3 * 1000;  d4 = d4 / 100;
                d5 = p - d1 * 100000 - d2 * 10000 - d3 * 1000 - d4 * 100; d5 = d5 / 10;
                d6 = p - d1 * 100000 - d2 * 10000 - d3 * 1000 - d4 * 100  - d5 * 10;

                // 612345
                // 561234
                // 456123
                // 345612
                // 234561
                r  = 100000*d6 + 10000*d1 + 1000*d2 + 100*d3 + 10*d4 + d5;
                r1 = 100000*d5 + 10000*d6 + 1000*d1 + 100*d2 + 10*d3 + d4;
                r2 = 100000*d4 + 10000*d5 + 1000*d6 + 100*d1 + 10*d2 + d3;
                r3 = 100000*d3 + 10000*d4 + 1000*d5 + 100*d6 + 10*d1 + d2;
                r4 = 100000*d2 + 10000*d3 + 1000*d4 + 100*d5 + 10*d6 + d1;
                if (is_prime(r))
                if (is_prime(r1))
                if (is_prime(r2))
                if (is_prime(r3))
                if (is_prime(r4))
                {
                    vmap[p] = true;
                    vmap[r] = true;
                    vmap[r1] = true;
                    vmap[r2] = true;
                    vmap[r3] = true;
                    vmap[r4] = true;

                    std::cout << p << std::endl;
                }
            }

        }
    }
    return vmap.size();
}


long long fact(long long n)
{
    long long r = 1;
    for(long long p=1;p<=n;p++) r = r * p;
    return r;
}

long long Euler034(long long N)
{
    long long sum = 0;
    long long sumf = 0;
    long long t;
    long long digit;
    for(long long p=3;p<=N;p++)
    {
        t = p;
        sumf = 0;
        while (t > 0)
        {
            digit = t % 10;
            t = t - digit; t = t/10;
            sumf+=fact(digit);
        }
        if (sumf == p)
        {
            sum+=p;
            std::cout << p << "  " << sumf << std::endl;
        }
    }
    return sum;
}

long long Euler033(long long N)
{
    long long digit1;
    long long digit2;
    std::map<long long , long long> vmap;
    for (long long a = 10; a<N; a++)
    for (long long b = a+1; b<N; b++)
    {
        // common
        long long aa = a;
        long long bb = b;
        long long p = 10;
        {
            while ((aa%p==0) && (bb%p==0))
            {
                aa = aa/p; bb = bb/p;
            }
        }

        if ((aa > 9) && (bb > 9))
        {
            digit1 = aa / 10;
            digit2 = bb % 10;
            if (digit1 == digit2)
            {
                long long anew = aa-10*digit1;
                long long bnew = (bb - digit2)/10;
                if (a*bnew == b*anew)
                {
                    vmap[a]=b;
                    std::cout << aa << " * " << bb << " " << anew << " "<< bnew << " "<< std::endl;
                }
            }

            digit1 = aa % 10;
            digit2 = bb / 10;
            if (digit1 == digit2)
            {
                long long anew = (aa - digit1)/10;
                long long bnew = bb - 10*digit2;
                if (a*bnew == b*anew)
                {
                    vmap[a]=b;
                    std::cout << aa << " * " << bb << " " << anew << " "<< bnew << " "<< std::endl;
                }
            }
        }
    }
    long long prda=1;
    long long prdb=1;
    for (auto & [aa,bb] : vmap)
    {
        prda *= aa;
        prdb *= bb;
    }
    for(long long p=2;p<N/2;p++)
    {
        while ((prda%p==0) && (prdb%p==0))
        {
            prda = prda/p; prdb = prdb/p;
        }
    }

    return prdb;
}

long long Euler032(long long N)
{
    //45228
    std::map<long long , bool> vmap;
    long long c;
    long long a;
    long long b;
    long long t;
    long long digit;
    bool ok;
    std::map<long long, long long> vdigit;
    for (a = 1; a<=N; a++)
    for (b = a+1; b<=N; b++)
    {
        ok = true;
        vdigit.clear();

        t = a;
        while (t > 0)
        {
            digit = t % 10;
            t = t - digit; t = t/10;
            if (vdigit.find(digit) != vdigit.end())
            {
                ok = false;
                break;
            };
            if (digit==0)
            {
                ok = false;
                break;
            };
            vdigit[digit] = 1;
        }

        if (ok)
        {
            t = b;
            while (t > 0)
            {
                digit = t % 10;
                t = t - digit; t = t/10;
                if (vdigit.find(digit) != vdigit.end())
                {
                    ok = false;
                    break;
                };
                if (digit==0)
                {
                    ok = false;
                    break;
                };
                vdigit[digit] = 1;
            }
        }

        if (ok)
        {
            t = a*b;
            c = t;
            while (t > 0)
            {
                digit = t % 10;
                t = t - digit; t = t/10;
                if (vdigit.find(digit) != vdigit.end())
                {
                    ok = false;
                    break;
                };
                if (digit==0)
                {
                    ok = false;
                    break;
                };
                vdigit[digit] = 1;
            }

            if (ok)
            if (vdigit.size() == 9)
            {
                vmap[c]=true;
                std::cout << a << " " << b << " " << a*b << std::endl;
            }
        }
    }
    long long sum=0;
    for (auto & [aa,bb] : vmap) sum += aa;
    return sum;
}

long long Euler031(long long N)
{
    //73682
    long long target  = N;
    long long ways = 0;

    for (long long a = target; a >= 0; a -= 200) {
        for (long long b = a; b >= 0; b -= 100) {
            for (long long c = b; c >= 0; c -= 50) {
                for (long long d = c; d >= 0; d -= 20) {
                    for (long long e = d; e >= 0; e -= 10) {
                        for (long long f = e; f >= 0; f -= 5) {
                            for (long long g = f; g >= 0; g -= 2) {
                                ways++;
                            }
                        }
                    }
                }
            }
        }
    }
    return ways;
}

long long Euler030(long long N)
{
    //443839
    std::map<uinteger_t, bool> vmap;
    uinteger_t p;
    uinteger_t sum;
    long long digit;
    long long t;

    for (long long i = 2; i <= N; i++)
    {
        sum = 0;
        t = i;
        while (t > 0)
        {
            digit = t % 10;
            t = t - digit;
            t = t/10;
            p = upow(digit,5);
            sum += p;
            if (sum > i) break;
            if (p > i) break;
        }
        if (sum == i) vmap[i]=true;
    }

    sum = 0;
    for(auto& [n, b] : vmap)
        sum += n;
    return (long long)sum;
}

long long Euler029(long long aa, long long bb)
{
    std::map<uinteger_t, bool> vmap;
    uinteger_t p;

    for (long long a = 2; a <= aa; a++)
    {
        for (long long b = 2; b <= bb; b++)
        {
            p = upow(a,b);
            vmap[p]=true;
        }
    }
    return vmap.size();
}

long long sum_spiral_corner(long long n)
{
    if (n==1) return 1;
    long long delta = n-1;
    return 4*n*n - 6*delta;
}

long long Euler028(long long N)
{
    //669171001
    long long sum = 0;
    for (long long i = 1; i <= N; i+=2)
    {
        sum+=sum_spiral_corner(i);
    }
    return sum;
}

long long number_primes_quadratic(long long a, long long b)
{
    long long cnt = 0;long long t;long long i = 0;
    while (true)
    {
        t = i*i + a*i + b;
        if (t<2) break;
        if (is_prime(t) == false)
            break;
        i++;
        cnt++;
    }
    return cnt;
}
long long Euler027(long long N)
{
    long long maxv = 0;long long maxa = 0;long long maxb = 0;long long cnt;
    for (long long a = -1*N + 1; a < N; a++)
    {
        for (long long b = -1*N; b <= N; b++)
        {
            cnt = number_primes_quadratic(a, b);
            if (cnt > maxv)
            {
              maxv = cnt;
              maxa = a;
              maxb = b;
              std::cout << a << " " << b<< " " << cnt << std::endl;
            }
        }
    }
    return maxa*maxb;
}


long long recu(long long  i, long long v, long long cnt, long long N)
{
    if (v == 1) return cnt;
    if (cnt > N) return cnt;
    v *= 10;
    v %= i;
    cnt++;
    return recu(i,v, cnt, N);
}

long long Euler026(long long N)
{
    long long maxv = 0;
    long long max_idx = 0;
    for (long long d = 2; d < N; d++)
    {
        long long cnt = 0;
        long v = 10%d;

        cnt = recu(d, v, cnt, N);
        if (cnt > maxv && cnt < N)
        {
          max_idx = d;
          maxv = cnt;
        }
    }
    return max_idx;
}

long long Euler025(size_t N)
{
    //4782
    long long idx = 3;
    uinteger_t fibo1 = 1;
    uinteger_t fibo2 = 1;
    uinteger_t t;
    uinteger_t lim = 1;
    for (size_t i = 0; i < N-1; i++) lim *= 10;
    while(true)
    {
        t = fibo1 + fibo2;
        fibo1 = fibo2;
        fibo2 = t;
        if (t >= lim) return idx;
        idx++;
    }
    return idx;
}

long long Euler024(long N)
{
    //2783915460
    std::vector<int> fac = {1, 1, 2, 6, 24, 120, 720, 5040, 40320, 362880};
    std::vector<int> arr = {0, 1, 2, 3, 4,  5,   6,   7,    8,     9};
    long long ret = 0;

    N--;
    for (long long  i = arr.size() - 1; i >= 0; i--)
    {
        int t = N / fac[i];
        N %= fac[i];
        ret = 10 * ret + arr[t];
        std::cout << i  << " " << t << " " << N << " " << ret << std::endl;
        arr.erase(arr.begin()+t);
    }
    std::cout << "return " << ret << std::endl;
    return ret;
}




const unsigned int AbundantLIMIT = 28124;
std::set<unsigned int> abundant;

unsigned int sumDivisor(unsigned int x)
{
    unsigned int sum = 1;
    unsigned int opposite;
    for (unsigned int div = 2; div * div <= x; div++)
    {
        if (x % div == 0)
        {
          sum += div;
          opposite = x / div;
          if (opposite != div) sum += opposite;
        }
    }
    return sum;
}

bool isAbundant(unsigned int x)
{
    if (x >= AbundantLIMIT) return true;
    for (auto i : abundant)
    {
        if (i >= x) return false;
        if (abundant.count(x - i) == 0) continue;
        return true;
    }
    return false;
}

long long Euler023()
{
    //4179871
    for (unsigned int i = 1; i < AbundantLIMIT; i++)
    {
        if (sumDivisor(i) > i)
          abundant.insert(i);
    }

    unsigned long long sum = 0;
    for (unsigned int i = 0; i < AbundantLIMIT; i++)
    {
        if (!isAbundant(i))
          sum += i;
    }
    return sum;
}



long long Euler022()
{
    long long sum = 0;
    std::ifstream is("names.txt", std::ifstream::in);
    if (is)
    {
        // get length of file:
        is.seekg (0, is.end);
        int length = (int)is.tellg();
        is.seekg (0, is.beg);

        char * buffer = new char [length];

        //std::cout << "Reading " << length << " characters... "<< std::endl;;
        is.read (buffer,length);
        std::string s;
        std::vector<std::string> v;

        if (is)
        {
            //"MARY","PATRICIA","LINDA","BARBARA","ELI
            for (int j=0; j < length; j++)
            {
                if (buffer[j]=='\"')
                {
                    if (s.size()>0)
                    {
                        v.push_back(s);
                        s.clear();
                    }
                }
                else if (buffer[j]!=',') s=s+buffer[j];
            }
            std::cout  << v.size() << std::endl;
        }
        is.close();

        std::sort(v.begin(), v.end());
        int sc = 0;
        for (size_t i=0; i < v.size(); i++)
        {
            s = v[i];
            sc = 0;
            for (size_t j=0; j < s.size(); j++)
            {
                sc += (1 + s[j] - 'A');
            }
            sum += (i+1) * sc;
        }
    }
    return sum;
}


long long sum_proper_divisors(long long N)
{
    long long sum = 0;
    for (int i=1; i <= N/2; i++)
    {
        if (N % i == 0) sum+=i;
    }
    return sum;
}
long long Euler021(long long N)
{
    long long sum = 0;
    long long t;

    std::vector<long long> vdivisors(N+1,0);

    for (int i=0; i < N; i++)
    {
        t = sum_proper_divisors(i+1);
        vdivisors[i] = t;
    }

    std::map<long long, bool> map_div;
    for (int j=0; j < N; j++)
    {
        t = vdivisors[j];
        if ((t <= N) && (t>1) && (t!=j+1))
        {
            if (vdivisors[t-1] == j+1)
            {
                map_div[j] = true;
                map_div[t-1] = true;
            }
        }
    }
    for (auto& [a , b] : map_div)
    {
        sum += a+1;
    }
    return sum;
}


long long Euler020(long long N)
{
    uinteger_t n = 1;
    uinteger_t sum = 0;
    for (int i=1; i < N; i++)
        n *= i;

    while (n >= 10)
    {
        sum += n % 10;
        n = n/10;
    }
    sum += n;
    return (long long)sum;
}

int day_of_week(int y, int m, int d)
{
  static int t[] = {0, 3, 2, 5, 0, 3, 5, 1, 4, 6, 2, 4};
  y -= m < 3;
  return (y + y/4 - y/100 + y/400 + t[m-1] + d) % 7;
}

long long Euler019()
{
    int d = 1;
    long long sum = 0;

    for(int y=1901;y<=2000;y++)
    for (int m=1; m <= 12; m++)
    {
        int day = day_of_week(y, m, d);
        if (day == 0) sum++;
    }
    return sum;
}

void visit_triangle(long long N, long long i, long long j, mat::matrix& mmaxpaths, mat::matrix& mnumber, mat::matrix& mvisit)
{
    long long t1 = 0;
    long long t2 = 0;
    if (mvisit(i,j) == 0)
    {
        if (mnumber(i,j) == 0) {mvisit(i,j) = 1; return;}

        mvisit(i,j) = 1;
        if (i+1<=N-1)                if (mvisit(i+1,j  ) == 0) visit_triangle(N, i+1, j,   mmaxpaths, mnumber, mvisit);
        if ((i+1<=N-1) && (j+1<=N))  if (mvisit(i+1,j+1) == 0) visit_triangle(N, i+1, j+1, mmaxpaths, mnumber, mvisit);

        t1 = 0; if (i+1<=N-1)                  t1 = mmaxpaths(i+1,j );
        t2 = 0; if ((i+1<=N-1) && (j+1<=N-1))  t2 = mmaxpaths(i+1,j+1);
        mmaxpaths(i, j) = mnumber(i, j)  + std::max(t1,t2);
    }
}

long long Euler018()
{
    mat::matrix mnumber(15, 15, 0);
    mat::matrix mmaxpaths(15, 15, 0);
    mat::matrix mvisit(15, 15, 0);
//75
//95 64
//17 47 82
//18 35 87 10
//20 04 82 47 65
//19 01 23 75 03 34
//88 02 77 73 07 63 67
//99 65 04 28 06 16 70 92
//41 41 26 56 83 40 80 70 33
//41 48 72 33 47 32 37 16 94 29
//53 71 44 65 25 43 91 52 97 51 14
//70 11 33 28 77 73 17 78 39 68 17 57
//91 71 52 38 17 14 91 43 58 50 27 29 48
//63 66 04 68 89 53 67 30 73 16 69 87 40 31
//04 62 98 27 23 09 70 98 73 93 38 53 60 04 23

    uinteger_t sum = 0;
    std::ifstream is("18.txt", std::ifstream::in);
    if (is)
    {
        // get length of file:
        is.seekg (0, is.end);
        int length = (int)is.tellg();
        is.seekg (0, is.beg);

        char * buffer = new char [length];

        //std::cout << "Reading " << length << " characters... "<< std::endl;;
        // read data as a block:
        is.read (buffer,length);
        if (is)
        {
            long num;
            int pos = 0;
            for (int i=0; i < 15; i++)
            {
                num=0;
                for (int j=0; j < i+1; j++)
                {
                    num = 10*(buffer[pos] - '0')+ (buffer[pos+1] - '0');
                    pos+=3;
                    mnumber(i, j) = num;
                    std::cout << i << " " << j << " " << num << std::endl;
                }
            }
        }
        is.close();
    }
    visit_triangle(15, 0, 0, mmaxpaths, mnumber, mvisit);
    return mmaxpaths(0, 0);
}


std::string number_to_str(long long N)
{
    if (N==0) return "";
    if (N==1) return "one";
    if (N==2) return "two";
    if (N==3) return "three";
    if (N==4) return "four";
    if (N==5) return "five";
    if (N==6) return "six";
    if (N==7) return "seven";
    if (N==8) return "eight";
    if (N==9) return "nine";
    if (N==10) return "ten";
    if (N==11) return "eleven";
    if (N==12) return "twelve";
    if (N==13) return "thirteen";
    if (N==14) return "fourteen";
    if (N==15) return "fifteen";
    if (N==16) return "sixteen";
    if (N==17) return "seventeen";
    if (N==18) return "eighteen";
    if (N==19) return "nineteen";
    if (N==20) return "twenty";
    if (N==30) return "thirty";
    if (N==40) return "forty";
    if (N==50) return "fifty";
    if (N==60) return "sixty";
    if (N==70) return "seventy";
    if (N==80) return "eighty";
    if (N==90) return "ninety";

    if ((N>20) && (N<30)) return number_to_str(20) + "" + number_to_str(N-20);
    if ((N>30) && (N<40)) return number_to_str(30) + "" + number_to_str(N-30);
    if ((N>40) && (N<50)) return number_to_str(40) + "" + number_to_str(N-40);
    if ((N>50) && (N<60)) return number_to_str(50) + "" + number_to_str(N-50);
    if ((N>60) && (N<70)) return number_to_str(60) + "" + number_to_str(N-60);
    if ((N>70) && (N<80)) return number_to_str(70) + "" + number_to_str(N-70);
    if ((N>80) && (N<90)) return number_to_str(80) + "" + number_to_str(N-80);
    if ((N>90) && (N<100)) return number_to_str(90) + "" + number_to_str(N-90);

    if ((N>=100) && (N<1000))
    {
        long long h = N / 100;
        long long d = (N-h*100);
        if (N==100) return "onehundred";
        if (d==0) return number_to_str(h) + "hundred";
        return number_to_str(h) + "hundredand" + number_to_str(d);
    }
    if (N==1000) return "onethousand";
    return "";
}
long long Euler017(long long N)
{
    long long sum = 0;
    for (long long i=1; i <= N; i++)
    {
        sum += number_to_str(i).size();
    }
    return sum;
}

long long Euler016(long long N)
{
    uinteger_t t = 1;
    uinteger_t sum = 0;
    for (long long i=1; i <= N; i++) t*=2;

    uinteger_t digit;
    while ( t >= 10)
    {
        digit = t % 10;
        t = t / 10;
        sum += digit;
    }
    sum += t;
    return (long long)sum;
}

void visit(long long N, long long i, long long j, mat::matrix& m, mat::matrix& mvisit)
{
    long long t1 = 0;
    long long t2 = 0;
    if (mvisit(i,j) == 0)
    {
        mvisit(i,j) = 1;
        if (i-1>=0) if (mvisit(i-1,j) == 0) visit(N, i-1, j, m, mvisit);
        if (j-1>=0) if (mvisit(i,j-1) == 0) visit(N, i, j-1, m, mvisit);
        t1 = 0; if (i-1>=0) t1 = m(i-1, j);
        t2 = 0; if (j-1>=0) t2 = m(i, j-1);
        m(i, j) = t1 + t2;
    }
}

long long Euler015(long long N)
{
    mat::matrix mpaths(N+1, N+1, 0);
    mat::matrix mvisit(N+1, N+1, 0);

    mvisit(0, 0) = 1;
    mpaths(1, 0) = 1; mvisit(1, 0) = 1;
    mpaths(0, 1) = 1; mvisit(0, 1) = 1;
    mpaths(1, 1) = 2; mvisit(1, 1) = 1;

    visit(N+1, N, N, mpaths, mvisit);
    return mpaths(N,N);
}


long long Euler014(long long N)
{
    std::map<long long , long long > map_sequ;
    long long max_len = 0;
    long long max_idx = 0;
    long long  cnt;
    for (int i = 1; i < N; i++)
    {
        long long s= i;
        cnt = 1;
        while(s!=1)
        {
            if (s%2==0) s = s/2;
            else s = 3*s+1;
            cnt++;

            if (map_sequ.find(s) != map_sequ.end())
            {
                cnt += map_sequ[s];
                cnt -= 1;
                break;
            }
        }
        if (map_sequ.find(i) == map_sequ.end())
            map_sequ[i] = cnt;

        if (cnt > max_len) {max_len = cnt; max_idx = i;}
        std::cout << i << " " << cnt << std::endl;
    }
    return max_idx;
}

long long Euler013()
{
    uinteger_t sum = 0;
    std::ifstream is("big100.txt", std::ifstream::in);
    if (is)
    {
        // get length of file:
        is.seekg (0, is.end);
        int length = (int)is.tellg();
        is.seekg (0, is.beg);

        char * buffer = new char [length];

        //std::cout << "Reading " << length << " characters... "<< std::endl;;
        // read data as a block:
        is.read (buffer,length);
        if (is)
        {
            uinteger_t num;
            for (int i=0; i < 100; i++)
            {
                // 50 digit;
                num=0;
                for (int j=0; j < 50; j++)
                {
                    num = 10*num + (buffer[51*i+j] - '0');
                }
                std::cout  << i << " " << num << std::endl;
                sum+=num;
            }
        }
        is.close();
    }
    uinteger_t m = 100*1000; m *= 100*1000 ;
    uinteger_t t =  sum;
    while( t > m) t = t/10;
    return (long long)t;
}

long long number_divisor(long long n)
{
    long long cnt = 0;
    for(long long i=1;i<=n;++i)
    {
        if(n%i==0) cnt++;
    }
    return cnt;
}


long long triangle_number(long long n)
{
    return n*(n+1)/2;
}
long long Euler012(long long N)
{
    for(long long i=1;i<=1000*1000;++i)
    {
        long long t = triangle_number(i);
        if (countDivisors(t) > N)
            return t;
    }
    return 0;
}

long long Euler011()
{
    std::string s =
    "08 02 22 97 38 15 00 40 00 75 04 05 07 78 52 12 50 77 91 08 "
    "49 49 99 40 17 81 18 57 60 87 17 40 98 43 69 48 04 56 62 00 "
    "81 49 31 73 55 79 14 29 93 71 40 67 53 88 30 03 49 13 36 65 "
    "52 70 95 23 04 60 11 42 69 24 68 56 01 32 56 71 37 02 36 91 "
    "22 31 16 71 51 67 63 89 41 92 36 54 22 40 40 28 66 33 13 80 "
    "24 47 32 60 99 03 45 02 44 75 33 53 78 36 84 20 35 17 12 50 "
    "32 98 81 28 64 23 67 10 26 38 40 67 59 54 70 66 18 38 64 70 "
    "67 26 20 68 02 62 12 20 95 63 94 39 63 08 40 91 66 49 94 21 "
    "24 55 58 05 66 73 99 26 97 17 78 78 96 83 14 88 34 89 63 72 "
    "21 36 23 09 75 00 76 44 20 45 35 14 00 61 33 97 34 31 33 95 "
    "78 17 53 28 22 75 31 67 15 94 03 80 04 62 16 14 09 53 56 92 "
    "16 39 05 42 96 35 31 47 55 58 88 24 00 17 54 24 36 29 85 57 "
    "86 56 00 48 35 71 89 07 05 44 44 37 44 60 21 58 51 54 17 58 "
    "19 80 81 68 05 94 47 69 28 73 92 13 86 52 17 77 04 89 55 40 "
    "04 52 08 83 97 35 99 16 07 97 57 32 16 26 26 79 33 27 98 66 "
    "88 36 68 87 57 62 20 72 03 46 33 67 46 55 12 32 63 93 53 69 "
    "04 42 16 73 38 25 39 11 24 94 72 18 08 46 29 32 40 62 76 36 "
    "20 69 36 41 72 30 23 88 34 62 99 69 82 67 59 85 74 04 36 16 "
    "20 73 35 29 78 31 90 01 74 31 49 71 48 86 81 16 23 57 05 54 "
    "01 70 54 71 83 51 54 69 16 92 33 48 61 43 52 01 89 19 67 48 ";

    long long prd = 0;
    long long maxprd = 0;
    std::vector<int> v;
    for(size_t i=0;3*i+1<s.size() ;i++)
    {
        v.push_back( 10*(s[3*i]-'0')  + 1*(s[3*i+1]-'0') );
    }

    for(size_t i=0;i<v.size() ;i++)
    {
        prd = 1; for(size_t j=0; (j<4) && (i+j < v.size());j++) prd*= v[i+j]; if (prd > maxprd) maxprd=prd; // right
        prd = 1; for(size_t j=0; (j<4) && (i + 20*j < v.size()) ;j++) prd*= v[i + 20*j]; if (prd > maxprd) maxprd=prd; // down
        prd = 1; for(size_t j=0; (j<4) && (i + j + 20*j < v.size()) ;j++) prd*= v[i + j + 20*j]; if (prd > maxprd) maxprd=prd; // diag
        prd = 1; for(size_t j=0; (j<4) && (i - j + 20*j < v.size())  ;j++) prd*= v[i - j + 20*j]; if (prd > maxprd) maxprd=prd; // diag

        if (prd > maxprd) maxprd=prd;
    }
    return maxprd;
}

long long Euler010(long N)
{
    long long sum = 0;
	for(long long i = 2; i < N  ; i++)
	{
		if(bitarray[i] == 0) // prime
		{
            sum+=i;

            if (i > last_index_processed)
            {
                for(long long j = i+1; j< SIZE_BITSET ; j++)
                {
                    if ( j % i == 0 ) // j a MUILTIPLE of prime i
                    {
                        bitarray[j] = 1; // MUILTIPLE of prime
                    }
                }
			}
		}
		last_index_processed = std::max(i, last_index_processed.load());
	}
	return sum;
}

long long Euler009()
{
    for(int i=1;i<1000 ;i++)
    {
        for(int j=i+1;j<1000 ;j++)
        {
            for(int k=j+1;k<1000 ;k++)
            {
                if ((i*i + j*j == k*k) && (i+j+k==1000))
                {
                    return i*j*k;
                }
            }
        }
    }
    return 0;
}

long long Euler008()
{
    std::string s =
    "73167176531330624919225119674426574742355349194934"
    "96983520312774506326239578318016984801869478851843"
    "85861560789112949495459501737958331952853208805511"
    "12540698747158523863050715693290963295227443043557"
    "66896648950445244523161731856403098711121722383113"
    "62229893423380308135336276614282806444486645238749"
    "30358907296290491560440772390713810515859307960866"
    "70172427121883998797908792274921901699720888093776"
    "65727333001053367881220235421809751254540594752243"
    "52584907711670556013604839586446706324415722155397"
    "53697817977846174064955149290862569321978468622482"
    "83972241375657056057490261407972968652414535100474"
    "82166370484403199890008895243450658541227588666881"
    "16427171479924442928230863465674813919123162824586"
    "17866458359124566529476545682848912883142607690042"
    "24219022671055626321111109370544217506941658960408"
    "07198403850962455444362981230987879927244284909188"
    "84580156166097919133875499200524063689912560717606"
    "05886116467109405077541002256983155200055935729725"
    "71636269561882670428252483600823257530420752963450";

    long long sum = 0;
    long long maxsum = 0;
    for(size_t i=13;i<s.size() ;i++)
    {
        sum = 1;
        for(size_t j=i-13;j<i ;j++)
        {
            sum *= (long) (s[j]-'0');
        }
        if (sum > maxsum) maxsum=sum;
    }
    return maxsum;
}


long long sum_square_diff(long num)
{
    long long  sum=0;
    long long sum_sq = 0;
    for(long long i = 1; i <= num ; i++)
    {
        sum+= i*i;
        sum_sq+=i;
    }
    long long sum_sqr = sum_sq * sum_sq;
    return sum_sqr-sum;
}

long long next_prime_factor(long long pp);
long long smallest_multiple(long long num)
{
    std::cout << "searching smallest_multiple 1 to " << num << std::endl;
    long long factor;
    long long p;
    std::map<long long , int> map_primes;
	for(long long i = 2; i <= num ; i++)
	{
        p = i;
        while (p > 1)
        {
            factor = next_prime_factor(p);
            if (factor > 1)
            {
                int cnt = 1;
                if (map_primes.find(factor) == map_primes.end())
                {
                    map_primes[factor] = 1;
                }
                else
                {
                    long long v = p/factor;
                    p = v;
                    if (v > 1)
                    {
                        while (next_prime_factor(v) == factor)
                        {
                            cnt++;
                            v = v/factor;
                            p = v;
                            if (factor == 1) break;
                        }
                    }
                    if (map_primes[factor] < cnt)
                    {
                        map_primes[factor] = cnt;
                    }
                }
            }
            else
            {
                p = 1;
            }
        }
    }

    p = 1;
    for(auto& [a, cnt] :  map_primes)
    {
        for(long long i = 0; i < cnt ; i++)
            p = p * a;
    }
    return p;
}

long long largest_palindrome(long long max_num)
{
    long long r = 0;
    long long mx= 0;
	for(long long i = 100; i < max_num ; i++)
	{
		for(long long j = 100; j < max_num ; j++)
        {
            r = i*j;
            if (is_palindrome(r))
            {
                if (r > mx) mx=r;
            }
        }
	}
	return mx;
}

long long find_prime(long long index_prime)
{
    std::cout << "searching prime index (Sieve of Eratosthenes algo) ... " << index_prime << std::endl;
    if (array_reset == false)
    {
        array_reset = true;
        bitarray.reset();
    }

	long long idx = 0;
	for(long long i = 2; i <= SIZE_BITSET  - 1 ; i++)
	{
        if (i % 1000 == 1)
            std::cout << "searching primes ... " << idx<< " : " << i << " SIZE_BITSET: " << SIZE_BITSET << std::endl;

		if(bitarray[i] == 0) // prime
		{
			idx++;
            if(idx == index_prime)
            {
                return i;
            }

            if (i > last_index_processed)
            {
                for(long long j = i+1; j< SIZE_BITSET ; j++)
                {
                    if ( j % i == 0 ) // j a MUILTIPLE of prime i
                    {
                        bitarray[j] = 1; // MUILTIPLE of prime
                    }
                }
			}
			last_index_processed = i;
		}
	}
	return 1;
}


long long largest_prime_factor(long long p)
{
    find_prime(10);

    long long last_factor = 1;
    long long pp = p;
    while (pp  > 1)
    {
        for(long long i = 2; i <= SIZE_BITSET - 1 ; i++)
        {
            if (i <= last_index_processed)
            {
                if (bitarray[i] == 0) // prime
                {
                    if (pp % i == 0)
                    {
                        last_factor = i;
                        pp = pp/ i;
                        if (pp == 1)
                            return last_factor;
                    }
                }
            }
            else
            {
                find_prime(last_index_processed+10);
            }
        }
    }
    return last_factor;
}

long long next_prime_factor(long long n)
{
    if (n == 1) return 1;
    if (n==2) return 2;

    if (array_reset == false)
    {
        array_reset = true;
        bitarray.reset();
    }

    long long next_factor = n;
    if (last_index_processed >= n)
    {
        if (n <= SIZE_BITSET)
        {
            if (bitarray[n] == 0) // prime
            {
                return n;
            }
        }
    }

    for(long long i = 2; i <= SIZE_BITSET - 1 ; i++)
    {
        if (i <= last_index_processed)
        {
            if (bitarray[i] == 0) // prime
            {
                if (n % i == 0)
                {
                    return i;
                }
            }
        }
        else
        {
            find_prime(last_index_processed+10);
        }
    }

    return next_factor;
}

long long mainEuler002()
{
    long long sum = 2;
    long long v1 = 1;
    long long v2 = 2;
    long long t;
    while(true)
    {
        t =  v1 + v2;
        v1 = v2;
        v2 = t;
        if (t%2 == 0) sum+=t;
        if (t >= 4e6) break;
    }
    return sum;
}

int mainEuler001()
{
    int r=0;
    for(int i=0;i<1000;i++)
    {
        if ((i%3==0) || (i%5==0)) r+=i;
    }
    return r;
}

int main()
{
    long long n = 0; std::string s;
    if (is_prime(2) != true) std::cout << "ERROR in is_prime(2) " << 2 << std::endl;
    if (is_prime(29) != true) std::cout << "ERROR in is_prime(29) " << 29 << std::endl;
    if (is_prime(29*17) != false) std::cout << "ERROR in is_prime(29*17) " << 29*17 << std::endl;
    if (is_prime(961748941) != true) {std::cout << "ERROR in is_prime(961748941) " << 961748941<< std::endl; return 0;}

    {
        RationalNumberTest rnst;
        std::pair<int, bool>  p = rnst.unit_tests();
        if (p.second != true)
        {
            std::cout << "ERROR in RationalNumberTest" << p.first<< std::endl;
            return 0;
        }
    }

    //n = mainEuler001() ; to_file("Euler001", n);
    //std::cout << "Euler001 " << n << std::endl;

    //n = mainEuler002() ; to_file("Euler002", n);
    //std::cout << "Euler002 " << n << std::endl;

    //find_prime(40);
    //if(bitarray[2] != 0)  std::cout << "ERROR in prime list" << std::endl;
    //if(bitarray[3] != 0)  std::cout << "ERROR in prime list" << std::endl;
    //if(bitarray[5] != 0)  std::cout << "ERROR in prime list" << std::endl;
    //if(bitarray[19] != 0)  std::cout << "ERROR in prime list" << std::endl;
    //if(bitarray[29] != 0)  std::cout << "ERROR in prime list" << std::endl;

    //if (next_prime_factor(9) != 3) std::cout << "ERROR in next_prime_factor" << 9 << std::endl;
    //if (next_prime_factor(8) != 2) std::cout << "ERROR in next_prime_factor" << 8 << std::endl;

    //n = largest_prime_factor(600851475143) ; to_file("Euler003", n);
    //std::cout << "Euler003 " << n << std::endl;

    //n = largest_palindrome(1000) ; to_file("Euler004", n);
    //std::cout << "Euler004 " << n << std::endl;

    //n = smallest_multiple(10); if (n!=2520) std::cout << "ERROR in smallest_multiple(10) != 2520 " << n << std::endl;
    //n = smallest_multiple(20) ; to_file("Euler005", n);
    //std::cout << "Euler005 " << n << std::endl;

    //n = sum_square_diff(10) ; if (n!=2640) std::cout << "ERROR in sum square 10 " << n << std::endl;
    //n = sum_square_diff(100); to_file("Euler006", n);
    //std::cout << "Euler006 " << n << std::endl;

    //n = find_prime(6); if (n!=13) std::cout << "ERROR in prime indexing " << n << std::endl;
    //n = find_prime(10001); to_file("Euler007", n);
    //std::cout << "Euler007 " << n << std::endl;

    //n = Euler008(); to_file("Euler008", n);
    //std::cout << "Euler008 " << n << std::endl;

    //n = Euler009(); to_file("Euler009", n);
    //std::cout << "Euler009 " << n << std::endl;

    //n = Euler010(10); if (n!=17) std::cout << "ERROR in prime sum " << n << std::endl;
    //n = Euler010(2000000); to_file("Euler010", n);
    //std::cout << "Euler010 " << n << std::endl;

    //n = Euler011(); to_file("Euler011", n);
    //std::cout << "Euler011 " << n << std::endl;

    //n = Euler012(500); to_file("Euler012", n);
    //std::cout << "Euler012 " << n << std::endl;

    //n = Euler013(); to_file("Euler013", n);
    //std::cout << "Euler013 " << n << std::endl;

    //n = Euler014(1000*1000); to_file("Euler014", n);
    //std::cout << "Euler014 " << n << std::endl;

    //n = Euler015(2); if (n!=6) std::cout << "ERROR in Euler015 " << n << std::endl;
    //n = Euler015(20); to_file("Euler015", n);
    //std::cout << "Euler015 " << n << std::endl;
//
//    n = Euler016(15); if (n!=26) std::cout << "ERROR in Euler016 " << n << std::endl;
//    n = Euler016(1000); to_file("Euler016", n);
//    std::cout << "Euler016 " << n << std::endl;
//
//    n = number_to_str(342).size(); if (n!=23) std::cout << "ERROR in number_to_str " << n << std::endl;
//    n = number_to_str(115).size(); if (n!=20) std::cout << "ERROR in number_to_str " << n << std::endl;
//    n = Euler017(1000); to_file("Euler017", n);
//    std::cout << "Euler017 " << n << std::endl;
//
//    n = Euler018(); to_file("Euler018", n);
//    std::cout << "Euler018 " << n << std::endl;
//
//    n = Euler019(); to_file("Euler019", n);
//    std::cout << "Euler019 " << n << std::endl;
//
//    n = Euler020(100); to_file("Euler020", n);
//    std::cout << "Euler020 " << n << std::endl;

//    n = sum_proper_divisors(220); if (n!=284) std::cout << "ERROR in sum_proper_divisors " << n << std::endl;
//    n = sum_proper_divisors(284); if (n!=220) std::cout << "ERROR in sum_proper_divisors " << n << std::endl;
//    n = Euler021(10000); to_file("Euler021", n);
//    std::cout << "Euler021 " << n << std::endl;
//
//    n = Euler022(); to_file("Euler022", n);
//    std::cout << "Euler022 " << n << std::endl;
//
//    n = Euler023(); to_file("Euler023", n);
//    std::cout << "Euler023 " << n << std::endl;
//
//    n = Euler024(1000*1000); to_file("Euler024", n);
//    std::cout << "Euler024 " << n << std::endl;
//
//    n = Euler025(1000); to_file("Euler025", n);
//    std::cout << "Euler025 " << n << std::endl;

//    n = Euler026(1000); to_file("Euler026", n);
//    std::cout << "Euler026 " << n << std::endl;
//
//    n = Euler027(1000); to_file("Euler027", n);
//    std::cout << "Euler027 " << n << std::endl;
//
//    n = Euler028(5); if (n!=101) std::cout << "ERROR in Euler028 " << n << std::endl;
//    n = Euler028(1001); to_file("Euler028", n);
//    std::cout << "Euler028 " << n << std::endl;
//
//    n = Euler029(5,5); if (n!=15) std::cout << "ERROR in Euler029 " << n << std::endl;
//    n = Euler029(100,100); to_file("Euler029", n);
//    std::cout << "Euler029 " << n << std::endl; //9183
//
//    n = Euler030(1000000); to_file("Euler030", n);
//    std::cout << "Euler030 " << n << std::endl;

//    n = Euler031(200); to_file("Euler031", n);
//    std::cout << "Euler031 " << n << std::endl;
//
//    n = Euler032(10000); to_file("Euler032", n);
//    std::cout << "Euler032 " << n << std::endl;
//
//    n = Euler033(100); to_file("Euler033", n);
//    std::cout << "Euler033 " << n << std::endl;
//
//    n = Euler034(100000); to_file("Euler034", n);
//    std::cout << "Euler034 " << n << std::endl;
//
//    n = Euler035(1000000); to_file("Euler035", n);
//    std::cout << "Euler035 " << n << std::endl;
//
//    if (is_palindrome_base2(585)!=true) std::cout << "ERROR in is_palindrome_base2 " << 585 << std::endl;
//    n = Euler036(1000000); to_file("Euler036", n);
//    std::cout << "Euler036 " << n << std::endl;
//
//    if (is_truncatable_prime(3797)!=true) std::cout << "ERROR in is_truncatable_prime " << 3797 << std::endl;
//    n = Euler037(1000000); to_file("Euler037", n);
//    std::cout << "Euler037 " << n << std::endl;
//
//    if (is_pandigital_9(192, n) == false) std::cout << "ERROR in is_pandigital_9 " << 192 << std::endl;
//    n = Euler038(10000); to_file("Euler038", n);
//    std::cout << "Euler038 " << n << std::endl;
//
//    n = Euler039(1000); to_file("Euler039", n);
//    std::cout << "Euler039 " << n << std::endl;
//
//    if (d_k(9) != 9) std::cout << "ERROR in d_k 9" << std::endl;
//    if (d_k(10) != 1) std::cout << "ERROR in d_k 10" << std::endl;
//    if (d_k(11) != 0) std::cout << "ERROR in d_k 11" << std::endl;
//    if (d_k(9+2*90 + 1) != 1) std::cout << "ERROR in d_k 9+2*90 + 1" << std::endl;
//    n = Euler040(); to_file("Euler040", n);
//    std::cout << "Euler040 " << n << std::endl;

//    n = Euler041(10*1000000); to_file("Euler041", n);
//    std::cout << "Euler041 " << n << std::endl;

//    n = Euler042(); to_file("Euler042", n);
//    std::cout << "Euler042 " << n << std::endl;
//
//    n = Euler043(); to_file("Euler043", n);
//    std::cout << "Euler043 " << n << std::endl;
//
//    n = Euler044(100000); to_file("Euler044", n);
//    std::cout << "Euler044 " << n << std::endl;

//    n = Euler045(100000); to_file("Euler045", n);
//    std::cout << "Euler045 " << n << std::endl;

//    n = Euler046(100000); to_file("Euler046", n);
//    std::cout << "Euler046 " << n << std::endl;

//    n = Euler047(1000000); to_file("Euler047", n);
//    std::cout << "Euler047 " << n << std::endl;
//
//    n = Euler048(1000); to_file("Euler048", n);
//    std::cout << "Euler048 " << n << std::endl;
//
//    n = Euler049(10000); to_file("Euler049", n);
//    std::cout << "Euler049 " << n << std::endl;
//
//    n = Euler050(1000000); to_file("Euler050", n);
//    std::cout << "Euler050 " << n << std::endl;
//
//    if (replace_digit(9334, 3, 5) != 9554) std::cout << "ERROR in replace_digit " << std::endl;
//    n = Euler051(1000000, 8); to_file("Euler051", n);
//    std::cout << "Euler051 " << n << std::endl;
//
//    n = Euler052(1000000); to_file("Euler052", n);
//    std::cout << "Euler052 " << n << std::endl;
//
//    n = Euler053(100); to_file("Euler053", n);
//    std::cout << "Euler053 " << n << std::endl;
////
////    //SEE Euler054.cpp
////
//    if (is_lych(10677) != true) std::cout << "ERROR in is_lych 10677" << std::endl;
//    if (is_lych(349) != false) std::cout << "ERROR in is_lych 349" << std::endl;
//    if (is_lych(196) != true) std::cout << "ERROR in is_lych 196" << std::endl;
//    if (is_lych(4994) != true) std::cout << "ERROR in is_lych 4994" << std::endl;
//    n = Euler055(10000); to_file("Euler055", n);
//    std::cout << "Euler055 " << n << std::endl;
//
//    n = Euler056(100); to_file("Euler056", n);
//    std::cout << "Euler056 " << n << std::endl;

    //n = Euler057(1000); to_file("Euler057", n);
    //std::cout << "Euler057 " << n << std::endl;

    //dec101_t d = Euler826(1000*1000); to_file("Euler826", d);
    //std::cout << "Euler826 " << d << std::endl;

    //n = Euler058(100000); to_file("Euler058", n);
    //std::cout << "Euler058 " << n << std::endl;

    //n = Euler059(); to_file("Euler059", n);
    //std::cout << "Euler059 " << n << std::endl;

//    if (concate(101, 54) != 10154) std::cout << "ERROR in concate" << std::endl;
//    n = Euler060(10000); to_file("Euler060", n);
//    std::cout << "Euler060 " << n << std::endl;

//    n = Euler061(10000); to_file("Euler061", n);
//    std::cout << "Euler061 " << n << std::endl;
//
//    n = Euler062(100000); to_file("Euler062", n);
//    std::cout << "Euler062 " << n << std::endl;
//
//    n = Euler063(1000); to_file("Euler063", n);
//    std::cout << "Euler063 " << n << std::endl;
//
//    n = Euler064(10000); to_file("Euler064", n);
//    std::cout << "Euler064 " << n << std::endl;
//
//    n = Euler065(100); to_file("Euler065", n);
//    std::cout << "Euler065 " << n << std::endl;
//
//    n = Euler066(1000); to_file("Euler066", n);
//    std::cout << "Euler066 " << n << std::endl;
//
//    n = Euler067(100); to_file("Euler067", n);
//    std::cout << "Euler067 " << n << std::endl;
//
//    std::string s = Euler068(5); to_file("Euler068", s);
//    std::cout << "Euler068 " << s << std::endl;
//
//    n = Euler069(1000000); to_file("Euler069", n);
//    std::cout << "Euler069 " << n << std::endl;
//
//    n = Euler070(10000000); to_file("Euler070", n);
//    std::cout << "Euler070 " << n << std::endl;

//    n = Euler071(1000000); to_file("Euler071", n);
//    std::cout << "Euler071 " << n << std::endl;

//    n = Euler072(1000000); to_file("Euler072", n);
//    std::cout << "Euler072 " << n << std::endl;

//    n = Euler073(12000); to_file("Euler073", n);
//    std::cout << "Euler073 " << n << std::endl;
//
//    n = Euler074(1000000); to_file("Euler074", n);
//    std::cout << "Euler074 " << n << std::endl;
//
//    n = Euler075(1500000); to_file("Euler075", n);
//    std::cout << "Euler075 " << n << std::endl;
//
//    n = Euler076(100); to_file("Euler076", n);
//    std::cout << "Euler076 " << n << std::endl;
//
//    n = Euler077(5000); to_file("Euler077", n);
//    std::cout << "Euler077 " << n << std::endl;
//
//    n = Euler078(60000); to_file("Euler078", n);
//    std::cout << "Euler078 " << n << std::endl;
//
//    n = Euler079(100000000); to_file("Euler079", n);
//    std::cout << "Euler079 " << n << std::endl;
//
//    n = Euler080(100); to_file("Euler080", n);
//    std::cout << "Euler080 " << n << std::endl;
//
//    n = Euler081(80); to_file("Euler081", n);
//    std::cout << "Euler081 " << n << std::endl;
//
//    s = Euler082(); to_file("Euler082", s);
//    std::cout << "Euler082 " << s << std::endl;
//
//    s = Euler083(); to_file("Euler083", s);
//    std::cout << "Euler083 " << s << std::endl;
//
//    // Euler084 Monopoly...
//    std::cout << "Euler084 TODO..." << std::endl;
//
//    s = Euler085(2000000); to_file("Euler085", s);
//    std::cout << "Euler085 " << s << std::endl;
//
//    n = Euler086(1000000); to_file("Euler086", n);
//    std::cout << "Euler086 " << n << std::endl;
//
//    n = Euler087(50000000); to_file("Euler087", n);
//    std::cout << "Euler087 " << n << std::endl;
//
//    s = Euler088(12000); to_file("Euler088", s);
//    std::cout << "Euler088 " << s << std::endl;
//
//    n = Euler089(1000); to_file("Euler089", n);
//    std::cout << "Euler089 " << n << std::endl;
//
//    n = Euler090(2); to_file("Euler090", n);
//    std::cout << "Euler090 " << n << std::endl;
//
//    n = Euler091(50); to_file("Euler091", n);
//    std::cout << "Euler091 " << n << std::endl;
//
//    n = Euler092(10000000); to_file("Euler092", n);
//    std::cout << "Euler092 " << n << std::endl;
//
//    n = Euler093(0); to_file("Euler093", n);
//    std::cout << "Euler093 " << n << std::endl;

//    n = Euler094(1000*1000*1000); to_file("Euler094", n);
//    std::cout << "Euler094 " << n << std::endl;
//
//    n = Euler095(1000*1000); to_file("Euler095", n);
//    std::cout << "Euler095 " << n << std::endl;

//    n = Euler096(50); to_file("Euler096", n);
//    std::cout << "Euler096 " << n << std::endl;

    // prime sieve
    //do_prime_sieve(30, true);
    prime_sieve_mt(30, true);

    n = Euler827(18); to_file("Euler827", n);
    std::cout << "Euler827 " << n << std::endl;

//    n = Euler828(); to_file("Euler828", n);
//    std::cout << "Euler828 " << n << std::endl;
//
//    n = Euler097(); to_file("Euler097", n);
//    std::cout << "Euler097 " << n << std::endl;
//
//    p098 c098;
//    n = c098.solve(); to_file("Euler098", n);
//    std::cout << "Euler098 " << n << std::endl;
//
//    p099 c099;
//    n = c099.solve(); to_file("Euler099", n);
//    std::cout << "Euler099 " << n << std::endl;
//
//    p100 c100;
//    n = c100.solve(); to_file("Euler100", n);
//    std::cout << "Euler100" << n << std::endl;


//void prime_sieve_mt(int max_number_of_thread)

    std::cout << "Done enter a number to exit " << std::endl;
    int a; std::cin >> a;
    return 0;
}
