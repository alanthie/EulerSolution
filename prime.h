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

    constexpr int SIZE_BITSET = (1 << 28) + 1; // > few million number scan (bits) < 64 MB
    std::bitset<SIZE_BITSET> bitarray;
    std::atomic<long long> last_index_processed = 0;
    bool array_reset = false;

    long long find_prime_in_bitarray(long long index_prime)
    {
        //std::cout << "searching prime index (Sieve of Eratosthenes algo) ... " << index_prime << std::endl;
        if (array_reset == false)
        {
            array_reset = true;
            bitarray.reset();
            std::cout << "Primes Eratosthenes array reset ... " << " SIZE_BITSET: " << SIZE_BITSET << std::endl;
        }

        long long idx = 0;
        for(long long i = 2; i <= SIZE_BITSET  - 1 ; i++)
        {
            if (i % 100 == 1)
                std::cout << "Primes Eratosthenes search ... " << "last_index_processed:" << last_index_processed << " SIZE_BITSET: " << SIZE_BITSET << std::endl;

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
            }
            last_index_processed = std::max(i, last_index_processed.load());
        }
        return 1;
    }

    //std::thread(fill_bitarray_primes))
    void fill_bitarray_primes()
    {
        long long index_prime = 0;
        while (last_index_processed < SIZE_BITSET - 1)
        {
            index_prime+=10000;
            find_prime_in_bitarray(index_prime);
        }
    }

    using prime_type = long long;
    using prime_container = std::vector<long long>;
    std::map<long long, bool> map_prime;
    std::map<uinteger_t, bool> map_uprime;

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

    bool is_uprime(uinteger_t n)
    {
        if (map_uprime.find(n) != map_uprime.end())
            return true;

        if (n<2) return false;
        if (n==2) {map_uprime[n]=true;return true;}
        if (n==666666667) {map_uprime[n]=true;return true;}

        if (n < last_index_processed.load())
        {
            long long i = (long long)n;
            bool b = (bitarray[i] == 0) ? true : false;
            return b;
        }

        uinteger_t N = n;
        for (uinteger_t i=2; i*i <= N; i++)
        {
            // a divisor exist
            if (N % i == 0) return false;
        }
        map_uprime[n]=true;
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

    uinteger_t next_uprime(uinteger_t p, uinteger_t n = 0)
    {
        if (n == 0)
        {
            uinteger_t i = p+1;
            while(true)
            {
                if (is_uprime(i)) return i;
                i++;
            }
            return 2;
        }

        uinteger_t N = n;
        for (uinteger_t i = p+1; i <= N; i++)
        {
            if (is_uprime(i)) return i;
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
        if (a < 0)
        {
            std::cout << "ERROR upow a < 0 " <<std::endl;
            return 1;
        }
        if ( b < 0)
        {
            std::cout << "ERROR upow b < 0 " <<std::endl;
            return 1;
        }
        uinteger_t r = 1;
        for (long long i = 1; i <= b; i+=1)
            r = r * a;
        return r;
    }

    uinteger_t power_modulo(uinteger_t a, uinteger_t power, uinteger_t mod)
    {
        // (a ⋅ b) mod m = [(a mod m) ⋅ (b mod m)] mod m
        if (power==0) return 1;
        if (power%2 == 1) return ((a % mod) * power_modulo(a, power-1, mod)) % mod;
        uinteger_t b = power_modulo(a, power/2, mod) % mod;
        return (b*b) % mod;
    }
    uinteger_t prod_modulo(uinteger_t a, uinteger_t b, uinteger_t mod)
    {
        // (a ⋅ b) mod m = [(a mod m) ⋅ (b mod m)] mod m
        uinteger_t r = ( a % mod) * ( b % mod);
        return r % mod;
    }
    uinteger_t sum_modulo(uinteger_t a, uinteger_t b, uinteger_t mod)
    {
        // (a ⋅ b) mod m = [(a mod m) ⋅ (b mod m)] mod m
        uinteger_t r = ( a % mod) + ( b % mod);
        return r % mod;
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

     std::vector<uinteger_t> uprime_factors(uinteger_t k)
     {
        std::vector<uinteger_t> r;
        uinteger_t t = k;
        for (uinteger_t i = 2; i <= t; i++)
        {
            if (t == 1) break;
            if (is_uprime(t))
            {
                r.push_back(t);
                t = 1;
            }
            else if (is_uprime(i))
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

    std::unordered_map<long long, long long> unique_prime_factors(long long k, size_t limit = 0)
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

                    if (limit > 0)
                        if (r.size() > limit) return r;
                }
            }
        }
        return r;
    }

    std::map<uinteger_t, long long> unique_uprime_factors(uinteger_t k, size_t limit = 0)
    {
        std::map<uinteger_t, long long> r;
        if (is_uprime(k)) {r[k]++; return r;}

        uinteger_t t = k;
        for (uinteger_t i = 2; i <= t; i++)
        {
            if (t == 1) break;
            if (is_uprime(i))
            {
                while(t % i == 0)
                {
                    //std::cout << "[prime :"<<i<<"]" << std::endl;
                    r[i]++;
                    t = t/i;

                    if (limit > 0)
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

    long long sumofFactors(long long n)
    {
        // Traversing through all prime factors.
        long long res = 1;
        long long k = n;
        for (long long i = 2; i <= sqrt(n); i++)
        {
            long long curr_sum = 1;
            long long curr_term = 1;
            while (n % i == 0)
            {
                // THE BELOW STATEMENT MAKES
                // IT BETTER THAN ABOVE METHOD
                //  AS WE REDUCE VALUE OF n.
                n = n / i;

                curr_term *= i;
                curr_sum += curr_term;
            }

            res *= curr_sum;
        }

        // This condition is to handle
        // the case when n is a prime
        // number greater than 2.
        if (n >= 2)
            res *= (1 + n);

        return res - k;
    }

    long long countDivisors(long long n)
    {
        long long cnt = 0;
        for (long long i = 1; i*i <= n; i++) {
            if (n % i == 0) {
                // If divisors are equal,
                // count only one
                if (n / i == i)
                    cnt++;

                else // Otherwise count both
                    cnt = cnt + 2;
            }
        }
        return cnt;
    }

    uinteger_t ucountDivisors(uinteger_t n)
    {
        uinteger_t cnt = 0;
        for (uinteger_t i = 1; i*i <= n; i++) {
            if (n % i == 0) {
                // If divisors are equal,
                // count only one
                if (n / i == i)
                    cnt++;

                else // Otherwise count both
                    cnt = cnt + 2;
            }
        }
        return cnt;
    }
}
#endif // PRIME_H_INCLUDED
