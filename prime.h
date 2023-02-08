#ifndef PRIME_H_INCLUDED
#define PRIME_H_INCLUDED

#include "common.h"
#include "matrix.hpp"

#include "atomic_bitvector.hpp"

namespace PRIME
{
    std::mutex mutex_prime_output;

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

    constexpr unsigned int SIZE_BITSET = (1 << 30) + 1; //  30 == scan 1 billion number (bits) < 128 MB

    // std::bitset<SIZE_BITSET> bitarray;
    // ATOMIC BITSET https://github.com/ekg/atomicbitvector
    atomicbitvector::atomic_bv_t bitarray(SIZE_BITSET);
    std::atomic<long long> last_index_processed = 0;

    bool        is_prime(long long n, long long istart);
    long long   next_prime(long long p, long long n, bool can_check_bitset_for_none_prime, long long istart);
    void        check_prime_sieve();
    bool        is_prime_brute_force(long long n);

    inline void fill_bitarray_primes_at_i(long long i)
    {
        if (i < SIZE_BITSET )
        {
            for(long long j = i*2; j< SIZE_BITSET ; j+=i)
            {
                bitarray.set(j); // MUILTIPLE of prime
            }
        }
    }

    class PrimeThread
    {
    public:

        std::atomic<bool> done = false;
        std::atomic<bool> has_work = false;
        std::atomic<bool> is_working = false;
        long long k;
        std::thread* t = nullptr;;

        PrimeThread() {}
        ~PrimeThread()
        {
            if (t!=nullptr)
            {
                t->join();

                delete t;
                t = nullptr;
            }
        }

        bool set_work(long long i)
        {
            if (has_work == false)
            {
                k = i;
                has_work = true;
                return true;
            }
            return false;
        }

        bool haswork() {return has_work.load();}

        void exit()
        {
            done = true;
        }

        void start(long long i)
        {
            if (t == nullptr)
            {
                t = new std::thread(&PrimeThread::run_loop, this);
                this->set_work(i);
            }
        }

        void run_loop()
        {
            while(done==false)
            {
                if (is_working == false)
                {
                    if (has_work == true)
                    {
                        is_working = true;
                        fill_bitarray_primes_at_i(k);

                        is_working = false;
                        has_work = false;
                    }
                    else
                    {
                        std::this_thread::sleep_for (std::chrono::microseconds(1));
                    }
                }
                else
                {
                    std::this_thread::sleep_for (std::chrono::microseconds(1));
                }
            }
        }
    };

    // COULD BE improve for speed
    std::vector<long long> get_n_next_prime(long long last_prime, size_t N, long long limit_prime)
    {
        std::vector<long long> r;
        long long entry = last_prime;
        long long lim =  1000+last_prime*last_prime;
        while(true)
        {
            last_prime = next_prime(last_prime, lim, true, (entry<2) ? 2 : entry);

            if (limit_prime >= 0)
                if (last_prime > limit_prime) break;

            if (r.size() >= N) break;
            r.push_back(last_prime);
        }
        return r;
    }

    void prime_sieve_mt(int max_number_of_thread, bool out)
    {
        std::vector<long long> vnext_primes;
        std::vector<std::thread> vt;
        long long last_prime = 1;

        auto tstart = std::chrono::steady_clock::now();

        while (true)
        {
            if (last_prime*last_prime >= SIZE_BITSET-1) break;

            vt.clear();
            vnext_primes = get_n_next_prime(last_prime, max_number_of_thread, SIZE_BITSET-1);

            if (vnext_primes.size() >= 1)
            {
                if (out)
                {
                    const std::lock_guard<std::mutex> lock(mutex_prime_output);
                    std::cout << "filling BITSET " <<  SIZE_BITSET << " with multiples of ";
                    for(size_t i = 0; i < vnext_primes.size() ; i++)
                    {
                        std::cout << vnext_primes[i] << " ";
                    }
                    std::cout << std::endl;
                }

                for(size_t i = 0; i < vnext_primes.size() ; i++)
                {
                    vt.push_back(std::thread(fill_bitarray_primes_at_i, vnext_primes[i]));
                }
                for (size_t k=0;k<vt.size();k++)  vt[k].join();

                last_index_processed = std::max(vnext_primes[vnext_primes.size() - 1], last_index_processed.load());
                last_prime = vnext_primes[vnext_primes.size() - 1];
            }
            else
            {
                break;
            }

            if (last_index_processed >= SIZE_BITSET - 1)
                break;
        }
        last_index_processed = SIZE_BITSET - 1;

        auto tend = std::chrono::steady_clock::now();
        //if (out)
        {
            const std::lock_guard<std::mutex> lock(mutex_prime_output);
            std::cout << "Elapsed time in milliseconds for prime sieve of BITSET size : "
                << SIZE_BITSET <<  " "
                << std::chrono::duration_cast<std::chrono::milliseconds>(tend - tstart).count()
                << " ms" << std::endl;
        }
        check_prime_sieve();
    }

    void prime_sieve_mt2(int max_number_of_thread, bool out)
    {
        bool first_time = true;
        std::vector<long long> vnext_primes;
        std::vector<PrimeThread*> vt;
        long long last_prime = 1;
        size_t cnt_done;
        double next_count(0);
        std::chrono::time_point<std::chrono::steady_clock> tstart_next ;
        std::chrono::time_point<std::chrono::steady_clock> tend_next ;

        auto tstart = std::chrono::steady_clock::now();

        while (true)
        {
            if (last_prime*last_prime >= SIZE_BITSET-1) break;

            tstart_next = std::chrono::steady_clock::now();
            vnext_primes = get_n_next_prime(last_prime, max_number_of_thread, SIZE_BITSET-1);
            tend_next  = std::chrono::steady_clock::now();
            next_count += (double)std::chrono::duration_cast<std::chrono::microseconds>(tend_next - tstart_next).count();

            if (vnext_primes.size() >= 1)
            {
                if (out)
                {
                    const std::lock_guard<std::mutex> lock(mutex_prime_output);
                    std::cout << "filling BITSET " <<  SIZE_BITSET << " with multiples of ";
                    for(size_t i = 0; i < vnext_primes.size() ; i++)
                    {
                        std::cout << vnext_primes[i] << " ";
                    }
                    std::cout << std::endl;
                }

                // First time
                if (vt.size() == 0)
                {
                    PrimeThread* t;
                    for(size_t i = 0; i < vnext_primes.size() ; i++)
                    {
                        t = new PrimeThread();
                        vt.push_back(t);
                        t->start(vnext_primes[i] );
                    }
                }

                if (first_time == false)
                {
                    for (size_t k=0;k<vnext_primes.size(); k++)
                    {
                        while (vt[k]->haswork() == true)
                        {
                            std::cout <<  "ERROR WAITING " <<  vnext_primes[k] << std::endl;
                            std::this_thread::sleep_for (std::chrono::microseconds(1));
                        }
                        vt[k]->set_work(vnext_primes[k]);
                    }
                }
                first_time = false;

                // wait all done
                while(true)
                {
                    cnt_done = 0;
                    for (size_t k=0;k<vnext_primes.size(); k++)
                    {
                        if (vt[k]->haswork() == false) cnt_done++;
                    }
                    if (cnt_done == vnext_primes.size()) break;

                    // sleep TODO
                    std::this_thread::sleep_for (std::chrono::microseconds(1));

                }

                last_index_processed = std::max(vnext_primes[vnext_primes.size() - 1], last_index_processed.load());
                last_prime = vnext_primes[vnext_primes.size() - 1];
            }
            else
            {
                break;
            }
            if (last_index_processed >= SIZE_BITSET - 1) break;
        }
        last_index_processed = SIZE_BITSET - 1;
        auto tend = std::chrono::steady_clock::now();

        for (size_t k=0;k<vt.size();k++)  vt[k]->exit();
        for (size_t k=0;k<vt.size();k++)  delete vt[k];

        //if (out)
        {
            const std::lock_guard<std::mutex> lock(mutex_prime_output);
            std::cout << "Elapsed time in milliseconds for prime sieve of BITSET size : "
                << SIZE_BITSET <<  " "
                << std::chrono::duration_cast<std::chrono::milliseconds>(tend - tstart).count()<< " ms"
                << " next total time " << next_count <<  " microseconds"
                << std::endl;
        }
        check_prime_sieve();
    }


    long long find_prime_in_bitarray(long long index_prime, bool out=false)
    {
        if (index_prime >= SIZE_BITSET - 1)
        {
            if (last_index_processed >= SIZE_BITSET - 1)
                return 1;
        }

        if (out)
        {
            const std::lock_guard<std::mutex> lock(mutex_prime_output);
            std::cout << "find_prime_in_bitarray " << index_prime << " SIZE_BITSET: " << SIZE_BITSET << " last_index_processed: " << last_index_processed<< std::endl;
        }


        long long idx = 0;
        std::vector<std::thread> vt;

        for(long long i = 2; i <= SIZE_BITSET - 1 ; i++)
        {
            if (out)
            if (i % 1000 == 0)
            {
                const std::lock_guard<std::mutex> lock(mutex_prime_output);
                std::cout << "Primes Eratosthenes search ... " << "i:" << i << " SIZE_BITSET: " << SIZE_BITSET << " last_index_processed: " << last_index_processed << std::endl;
            }

            if (bitarray.test(i) == false) // prime
            {
                idx++;
                if(idx == index_prime)
                {
                    return i;
                }

                if (i > last_index_processed)
                {
                    fill_bitarray_primes_at_i(i);
                }
            }
            last_index_processed = std::max(i, last_index_processed.load());
        }
        return 1;
    }

    //std::thread(fill_bitarray_primes))
    void fill_bitarray_primes(bool out=false)
    {
        long long index_prime = 0;
        while (last_index_processed < SIZE_BITSET - 1)
        {
            index_prime+=10000;
            find_prime_in_bitarray(index_prime, out);
        }
    }


    using prime_type = long long;
    using prime_container = std::vector<long long>;
    std::map<long long, bool> map_prime;
    std::map<uinteger_t, bool> map_uprime;

    bool is_prime(long long n, long long istart = 2)
    {
        if (map_prime.find(n) != map_prime.end())
            return true;

        if (n<2) return false;
        else if (n==2) return true;
        else if (n==9901) return true;
        else if (n==246209) return true;
        else if (n==7384709) return true;
        else if (n==27961) return true;
        else if (n==909091) return true;
        else if (n==5882353) return true;
        else if (n==215659) return true;
        else if (n==1000001) return false;
        else if (n==100000001) return false;
        else if (n==1000000001) return false;
        else if (n==66666667) return true;
        else if (n==666666667) return true;
        else if (n==6666666667) return false;
        else if (n==246209) return true;
        else if (n==7384709) return true;
        else if (n==5879) return true;
        else if (n==26801) return true;
        else if (n==39791) return true;
        else if (n==1873) return true;
        else if (n==41161) return true;
        else if (n==50867) return true;
        else if (n==985694468327) return true;

        if (n < last_index_processed.load())
        {
            bool b = (bitarray.test(n) == false) ? true : false;
            return b;
        }

        uinteger_t N = n;
        for (uinteger_t i=istart; i*i <= N; i++)
        {
            // a divisor exist
            if (N % i == 0) return false;
        }
        map_prime[n]=true;
        return true;
    }

    bool is_prime_brute_force(long long k)
    {
        if (k<2) return false;
        else if (k==2) return true;

        for (long long i=2; i*i <= k; i++)
        {
            // a divisor exist
            if (k % i == 0)
                return false;
        }
        return true;
    }

    bool is_uprime(uinteger_t n)
    {
        if (map_uprime.find(n) != map_uprime.end())
            return true;

        if (n<2) return false;
        else if (n==2) {map_uprime[n]=true;return true;}
        else if (n==666666667) {map_uprime[n]=true;return true;}

        if (n < last_index_processed.load())
        {
            long long i = (long long)n;
            bool b = (bitarray.test(i) == false) ? true : false;
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

    long long next_prime(long long p, long long n, bool can_check_bitset_for_none_prime = false, long long istart = 2)
    {
        long long N = n;
        for (long long i = p+1; i <= N; i++)
        {
            if (can_check_bitset_for_none_prime)
            {
                if (i <= SIZE_BITSET - 1)
                    if (bitarray.test(i) == true) // MUILTIPLE of prime
                        continue;
            }
            if (is_prime(i, istart))
                return i;
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

        //https://www.calculator.net/prime-factorization-calculator.html?cvar=20000000000001&x=40&y=17
        if (k==1818181818181) return {246209, 7384709};
        else if (k==181818181818181) return {29, 5879, 26801, 39791};
        else if (k==2000000000001) return {3, 43, 2347, 6605827};
        else if (k==20000000000001) return {3, 7, 5981, 159234401};
        else if (k==200000000000001) return {3, 17,  1873, 41161, 50867};
        else if (k==2000000000000001) return {3, 127, 138563, 37884167};
        else if (k==20000000000000001) return {3, 179, 12713, 2929595521};
        else if (k==200000000000000001) return {3, 44087, 691381,2187161};
        else if (k==285714285714285) return {3, 5, 952381,19999999 };
        else if (k==153846153846153) return {3, 3, 31, 197003, 662369};
        else if (k==9478672985781) return {3, 3, 719,32999, 44389};
        else if (k==44543429844097) return {823, 31663, 1709353};
        else if (k==2000000000000000001) return {3, 261382937, 2550536291};
        else if (k==19801980198019801) return {35837, 552556860173};
        else if (k==1980198019801) return {79, 12823, 1954753};

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
            if (is_prime(t))
            {
                r.push_back(t);
                t = 1;
            }
        }
        std::sort(r.begin(), r.end());
        return r;
     }

     std::vector<uinteger_t> uprime_factors(uinteger_t k, bool stop_ifbiggerprime = false, long long prime_limit = 0, std::vector<long long> vp = {})
     {
        std::vector<uinteger_t> r;
        uinteger_t t = k;
        for (uinteger_t i = 2; i <= t; i++)
        {
            if (t == 1) break;
            if (is_uprime(i))
            {
                while(t % i == 0)
                {
                    r.push_back(i);
                    t = t/i;;

                    if (stop_ifbiggerprime)
                    {
                        if (i >= prime_limit)
                        {
                            if (vp.size() == 3)
                                if ((i!=vp[0]) && (i!=vp[1]) && (i!=vp[2]))
                                    return r;

                            if (vp.size() == 4)
                                if ((i != vp[0]) && (i != vp[1]) && (i != vp[2]) && (i != vp[3]))
                                    return r;
                        }
                    }
                }

                if (is_uprime(t))
                {
                    r.push_back(t);
                    t = 1;
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

    std::vector<long long>  divisors(long long n)
    {
        std::vector<long long> v;
        if (n==2000002) return {1, 2, 101, 202, 9901, 19802, 1000001, 2000002};
        if (n==20000002) return {1, 2, 11, 22, 909091,  1818182,  100000001, 20000002};
        if (n==200000002) return {1, 2, 17, 34, 5882353, 11764706, 1000000001, 200000002};
        if (n==2000000002) return {1, 2, 7, 11, 13, 14, 19, 22, 26, 38, 77, 91, 133, 143, 154, 182, 209, 247, 266, 286, 418, 494, 1001, 1463, 1729, 2002, 2717, 2926, 3458, 5434, 19019, 38038, 52579, 105158, 368053, 578369, 683527, 736106, 999001, 1156738, 1367054, 1998002, 4048583, 4784689, 6993007, 7518797, 8097166, 9569378, 10989011, 12987013, 13986014, 15037594, 21978022, 25974026, 52631579, 76923077, 90909091, 105263158, 142857143, 153846154, 181818182, 285714286, 1000000001, 2000000002};
        if (n==20000000002) return {1, 2, 101, 202, 3541, 7082, 27961, 55922, 357641, 715282, 2824061, 5648122, 99009901, 198019802, 10000000001, 20000000002};
        if (n==200000000002) return {1, 2, 11, 22, 23, 46, 121, 242, 253, 506, 2783, 4093, 5566, 8186, 8779, 17558, 45023, 90046, 94139, 96569, 188278, 193138, 201917, 403834, 495253, 990506, 1035529, 1062259, 2071058, 2124518, 2221087, 4442174, 11390819, 22781638, 24431957, 35932447, 48863914, 71864894, 395256917, 790513834, 826446281, 1652892562, 4347826087, 8695652174, 9090909091, 18181818182, 100000000001, 200000000002};
        if (n==2000000000002) return {1, 2, 73, 137, 146, 274, 10001, 20002, 99990001, 199980002, 7299270073, 13698630137, 14598540146, 27397260274, 1000000000001, 2000000000002};
        if (n==20000000000002) return {1, 2, 11, 22, 859, 1718, 9449, 18898, 1058313049, 2116626098, 11641443539, 23282887078, 909090909091, 1818181818182, 10000000000001, 20000000000002};
        if (n==200000000000002) return {1, 2, 29, 58, 101, 202, 281, 562, 2929, 5858, 8149, 16298, 28381, 56762, 823049, 1646098, 121499449, 242998898, 3523484021, 7046968042, 12271444349, 24542888698, 34141345169, 68282690338, 355871886121, 711743772242, 990099009901, 1980198019802, 3448275862069, 6896551724138, 100000000000001, 200000000000002};
        if (n==2000000000000002) return {1, 2, 7, 11, 13, 14, 22, 26, 77, 91, 143, 154, 182, 211, 241, 286, 422, 482, 1001, 1477, 1687, 2002, 2161, 2321, 2651, 2743, 2954, 3133, 3374, 4322, 4642, 5302, 5486, 6266, 9091, 15127, 16247, 18182, 18557, 19201, 21931, 23771, 28093, 30173, 30254, 32494, 34463, 37114, 38402, 43862, 47542, 50851, 56186, 60346, 63637, 68926, 100001, 101702, 118183, 127274, 166397, 196651, 200002, 211211, 236366, 241241, 309023, 332794, 355957, 393302, 422422, 455971, 482482, 520801, 559361, 618046, 661063, 700007, 711914, 827281, 911942, 1041602, 1118722, 1300013, 1322126, 1400014, 1654562, 1918201, 2163161, 2190931, 2600026, 3191797, 3645607, 3836402, 3915527, 4326322, 4381862, 4627441, 5015681, 5728811, 5927623, 6383594, 6770413, 7271693, 7291214, 7831054, 9100091, 9254882, 10031362, 11457622, 11855246, 13427407, 13540826, 14543386, 15336517, 18200182, 19645651, 21100211, 24100241, 24936613, 26854814, 28482103, 30673034, 35109767, 39291302, 40101677, 41493361, 42200422, 47392891, 48200482, 49873226, 50901851, 56964206, 65203853, 70219534, 74474543, 80203354, 82986722, 94785782, 101803702, 109889011, 130407706, 137519557, 147701477, 148949086, 168701687, 174556291, 199374721, 216102161, 219778022, 255393463, 274302743, 275039114, 295402954, 313303133, 337403374, 349112582, 398749442, 432204322, 456426971, 462286441, 510786926, 521321801, 548605486, 626606266, 769223077, 912853942, 924572882, 1042643602, 1208779121, 1428557143, 1512715127, 1538446154, 1787754241, 1920119201, 2193121931, 2417558242, 2809328093, 2857114286, 3025430254, 3236005087, 3575508482, 3840238402, 4145232361, 4386243862, 4734601891, 5085150851, 5618656186, 6009723733, 6472010174, 8290464722, 8461453847, 9469203782, 9999900001, 10170301702, 12019447466, 15714128573, 16922907694, 19665296651, 19999800002, 29016626527, 31428257146, 33142213237, 35596055957, 39330593302, 42068066131, 45597555971, 52080620801, 53888020693, 58033253054, 61549824583, 66106961063, 66284426474, 71192111914, 84136132262, 91195111942, 104161241602, 107776041386, 109998900011, 123099649166, 132213922126, 219997800022, 319182891797, 364564345607, 377216144851, 430848772081, 462748727441, 592768227623, 638365783594, 677048070413, 729128691214, 754432289702, 861697544162, 925497454882, 999000999001, 1185536455246, 1354096140826, 1998001998002, 4149377593361, 4739336492891, 6993006993007, 8298755186722, 9478672985782, 10989010989011, 12987012987013, 13986013986014, 21978021978022, 25974025974026, 76923076923077, 90909090909091, 142857142857143, 153846153846154, 181818181818182, 285714285714286, 1000000000000001, 2000000000000002};
        if (n==20000000000000002) return {1, 2, 353, 449, 641, 706, 898, 1282, 1409, 2818, 69857, 139714, 158497, 226273, 287809, 316994, 452546, 497377, 575618, 632641, 903169, 994754, 1265282, 1806338, 24659521, 31365793, 44778337, 49319042, 62731586, 89556674, 98428513, 101596577, 196857026, 203193154, 223322273, 318818657, 405522881, 446644546, 637637314, 811045762, 11072124929, 15806752961, 20105473313, 22144249858, 31613505922, 34745265089, 40210946626, 44194402337, 63092676833, 69490530178, 88388804674, 126185353666, 143149576993, 286299153986, 7097232079489, 14194464158978, 15600624024961, 22271714922049, 28328611898017, 31201248049922, 44543429844098, 56657223796034, 10000000000000001, 20000000000000002};
        if (n==200000000000000002) return {1, 2, 11, 22, 103, 206, 1133, 2266, 4013, 8026, 44143, 88286, 413339, 826678, 4546729, 9093458, 21993833369, 43987666738, 241932167059, 483864334118, 2265364837007, 4530729674014, 24919013207077, 49838026414154, 88261253309797, 176522506619594, 970873786407767, 1941747572815534, 9090909090909091, 18181818181818182, 100000000000000001, 200000000000000002};
        if (n==2000000000000000002) return {1, 2, 101, 202, 9901, 19802, 1000001, 2000002, 999999000001, 1999998000002, 100999899000101, 201999798000202, 9900990099009901, 19801980198019802,  1000000000000000001,  2000000000000000002};
        //???
        //if (n==20000000000000000002) return {1, 2 , 11 , 909090909090909091, 10000000000000000001, 20000000000000000002};

        for (long long i = 1; i <= n/2; i++)
        {
            if (n % i == 0)
            {
                v.push_back(i);

//                // If divisors are equal,
//                // count only one
//                if (n / i == i)
//                    cnt++;
//
//                else // Otherwise count both
//                    cnt = cnt + 2;
            }
        }
        return v;
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


    void check_prime_sieve()
    {
        //auto tstart = std::chrono::steady_clock::now();
        //auto tend = std::chrono::steady_clock::now();
//        std::cout << "Elapsed time in milliseconds for prime sieve of BITSET size : "
//            << SIZE_BITSET <<  " "
//            << std::chrono::duration_cast<std::chrono::milliseconds>(tend - tstart).count()
//            << " ms" << std::endl;

        std::cout << "last_index_processed "<< last_index_processed << std::endl;
        std::cout << "SIZE_BITSET "<< SIZE_BITSET << std::endl;
        if (SIZE_BITSET >= 999983)
        {
            if(bitarray.test(999983) == false) std::cout << "ok 999983 is prime"<< std::endl;
            else std::cout << "ERROR 999983 should be prime"<< std::endl;
        }
        if (SIZE_BITSET >= 999983+1)
        {
            if(bitarray.test(999983+1) == true) std::cout << "ok 999983+1 is not prime"<< std::endl;
            else std::cout << "ERROR 999983+1 should not be prime"<< std::endl;
        }
        if (SIZE_BITSET >= 198491329)
        {
            if(bitarray.test(198491329) == false) std::cout << "ok 198491329 is prime"<< std::endl;
            else std::cout << "ERROR 198491329 should be prime"<< std::endl;
        }
        if (SIZE_BITSET >= 982451653)
        {
            if(bitarray.test(982451653) == false) std::cout << "ok 982451653 is prime"<< std::endl;
            else std::cout << "ERROR 982451653 should be prime"<< std::endl;
        }
        if (SIZE_BITSET >= 982451653+1)
        {
            if(bitarray.test(982451653+1) == true) std::cout << "ok 982451653+1 is not prime"<< std::endl;
            else std::cout << "ERROR 982451653+1 should not be prime"<< std::endl;
        }

        bool b; bool ba;
        for (long long  i = 2; i <= 10000000; i++)
        {
            if (i < last_index_processed.load())
            {
                b = is_prime_brute_force(i);
                if (bitarray.test(i) == false) ba = true; else ba = false;
                if (b!=ba) std::cout << "ERROR with " << i << " " << b << " " << ba << std::endl;
            }
        }
        std::cout << "ALL number to 10000000 checked " << std::endl;
    }
}
#endif // PRIME_H_INCLUDED
