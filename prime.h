#ifndef PRIME_H_INCLUDED
#define PRIME_H_INCLUDED

#include "common.h"
#include "matrix.hpp"

namespace PRIME
{
    class random_engine
    {
      public:
        std::random_device                      rd;
        std::mt19937                            mt;
        std::uniform_real_distribution<double>  dist;

        random_engine() : rd{}, mt{rd()}, dist{0.0, 1.0}
        {
            seed();
        }

        double get_rand()
        {
          return dist(mt);
        }

        void seed()
        {
            srand ((unsigned int)time(NULL));
            int n = rand() % 100;
            for (int i=0;i<n;i++) get_rand(); // random seed
        }
    };

    constexpr int SIZE_BITSET = (1 << 24) + 1; // > few million number scan (bits) < 64 MB
    std::bitset<SIZE_BITSET> bitarray;
    long long last_index_processed = 0;
    bool array_reset = false;

    using prime_type = long long;
    using prime_container = std::vector<long long>;
    std::map<long long, bool> map_prime;

    bool is_prime(long long n)
    {
        if (map_prime.find(n) != map_prime.end())
            return true;
        if (n<2) return false;
        if (n==2) return true;

        uinteger_t N = n;
        for (uinteger_t i=2; i*i <= N; i++)
        {
            // a divisor exist
            if (N % i == 0) return false;
        }
        map_prime[n]=true;
        return true;
    }

    long long next_prime(long long p, long long n)
    {
        long long N = n;
        for (long long i = p+1; i <= N; i++)
        {
            if (is_prime(i)) return i;
        }
        return 2;
    }

    prime_container list_prime(int N)
    {
        prime_container r;
        for (prime_type i=2; i <= N; i++)
        {
            if (is_prime(i))
                r.push_back(i);
        }
        return r;
    }

    void to_file(std::string title, long long n)
    {
        std::ofstream myfile;
        myfile.open ("solution.txt", std::ios_base::app);
        myfile << title  << ": " << n << '\n';
        myfile.close();
    }
    void to_file(std::string title, std::string n)
    {
        std::ofstream myfile;
        myfile.open ("solution.txt", std::ios_base::app);
        myfile << title  << ": " << n << '\n';
        myfile.close();
    }
    void to_file(std::string title, dec101_t n)
    {
        std::ofstream myfile;
        myfile.open ("solution.txt", std::ios_base::app);
        myfile << title  << ": " << std::setprecision(20) << n << '\n';
        myfile.close();
    }

    uinteger_t upow(long long a, long long b )
    {
        uinteger_t r = 1;
        for (long long i = 1; i <= b; i+=1)
            r = r * a;
        return r;
    }

    uinteger_t ufact(long long a )
    {
        uinteger_t r = 1;
        for (long long i = 2; i <= a; i+=1)
            r = r*i;
        return r;
    }

    uinteger_t ucombination(long long r, long long n )
    {
        uinteger_t a = ufact(r);
        uinteger_t b = ufact(n-r);
        uinteger_t c = ufact(n);
        return c / (a*b);
    }

    bool is_palindrome_base10(long long n)
    {
        int digit;
        long long r = n;
        long long reversal_number = 0;;
        while (r > 0)
        {
            digit = r % 10;
            reversal_number = 10 * reversal_number + digit;
            r = r/10;
        }
        return (reversal_number == n);
    }

    bool is_palindrome_base2(long long k)
    {
        std::vector<int> v;
        int digit;
        long long r = k;
        while (r > 0)
        {
            digit = r % 2;
            v.push_back(digit);
            r = r/2;
        }
        long long sz = v.size();
        for (long long i = 0; i < sz; i++)
        {
            if (v[i] != v[sz - i - 1]) return false;
        }
        return true;
    }

    bool is_truncatable_prime(long long k)
    {
        if (is_prime(k) == false) return false;

        std::vector<int> v;
        long long r = k;
        int digit;
        while (r > 0)
        {
            digit = r % 10;
            v.push_back(digit);
            r = r/10;
            if (r > 0)
                if (is_prime(r) == false) return false;
        }

        r = k;
        long long p;
        long long m;
        long long sz = v.size();
        for (long long i = sz; i > 0; i--)
        {
            p = (long long)upow(10, i-1);
            m = r / p;
            r = r - m*p;
            if (r > 0)
                if (is_prime(r) == false) return false;
        }
        return true;
    }

    long count_digit(int d, std::vector<int>& v)
    {
        long n = 0;
        for (size_t i = 0; i < v.size(); i++)
        {
            if (v[i] == d) n++;;
        }
        return n;
    }

    bool all_digit_unique(std::vector<int>& v)
    {
        for (int i = 0; i < 10; i++)
        {
            if (count_digit(i, v) > 1) return false;
        }
        return true;;
    }


    std::vector<int> digits10(long long k, bool keep_order = false)
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

        if (keep_order)
        {
            std::vector<int> vv;
            for (size_t i = 0; i < v.size(); i++)
            {
                vv.push_back(v[v.size() - i - 1]);
            }
            return vv;
        }
        return v;
    }

    std::vector<int> udigits10(uinteger_t k)
    {
        std::vector<int> v;
        uinteger_t r = k;
        uinteger_t digit;
        while (r > 0)
        {
            digit = r % 10;
            v.push_back((int)digit);
            r = r/10;
        }
        return v;
    }

    std::map<int, int> digitsmap10(long long k)
    {
        std::map<int, int> v;
        long long r = k;
        int digit;
        while (r > 0)
        {
            digit = r % 10;
            v[digit]++;
            r = r/10;
        }
        return v;
    }

    int digit_at(size_t pos, long long k)
    {
        // 01...n
        std::vector<int> v = digits10(k);
        if (pos < v.size())
            return v[v.size() - pos - 1];
        return 0;
    }

     std::vector<long long> prime_factors(long long k)
     {
        std::vector<long long> r;
        long long t = k;
        for (long long i = 2; i <= t; i++)
        {
            if (t == 1) break;
            if (is_prime(i))
            {
                while(t % i == 0)
                {
                    r.push_back(i);
                    t = t/i;;
                }
            }
        }
        return r;
     }

    std::unordered_map<long long, long long> unique_prime_factors(long long k, size_t limit)
    {
        std::unordered_map<long long, long long>r;
        if (is_prime(k)) {r[k]++; return r;}

        long long t = k;
        for (long long i = 2; i <= t; i++)
        {
            if (t == 1) break;
            if (is_prime(i))
            {
                while(t % i == 0)
                {
                    r[i]++;
                    t = t/i;

                    if (r.size() > limit) return r;
                }
            }
        }
        return r;
    }

    bool is_palindrome(long long n)
    {
        int digit;
        long long r = n;
        long long reversal_number = 0;;
        while (r > 0)
        {
            digit = r % 10;
            reversal_number = 10 * reversal_number + digit;
            r = r/10;
        }
        return (reversal_number == n);
    }

    bool is_upalindrome(uinteger_t n)
    {
        uinteger_t digit;
        uinteger_t r = n;
        uinteger_t reversal_number = 0;;
        while (r > 0)
        {
            digit = r % 10;
            reversal_number = 10 * reversal_number + digit;
            r = r/10;
        }
        return (reversal_number == n);
    }

    uinteger_t reverse_unumber(uinteger_t n)
    {
        uinteger_t digit;
        uinteger_t r = n;
        uinteger_t reversal_number = 0;
        while (r > 0)
        {
            digit = r % 10;
            reversal_number = 10 * reversal_number + digit;
            r = r/10;
        }
        return reversal_number;
    }

    dec101_t phi(long long n) {
        dec101_t result = n;
        for (long long i = 2; i * i <= n; i++) {
            if (n % i == 0) {
                while (n % i == 0)
                    n /= i;
                result -= result / i;
            }
        }
        if (n > 1)
            result -= result / n;
        return result;
    }

    long long longfact(long long n)
    {
        long long r = 1;
        for (long long p = 1; p <= n; p++) r = r * p;
        return r;
    }

    long long GCD(long long  a, long long  b)
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
    long long gcd(long long a, long long b)
    {
        long long c;
        while (b) {
            c = b;
            b = a % b;
            a = c;
        }
        return a;
    }

    long long to_long(std::string s)
    {
        return atoi(s.data());
    }
}

#endif // PRIME_H_INCLUDED
