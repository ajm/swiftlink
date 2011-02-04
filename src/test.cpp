using namespace std;

#include <cstdio>
#include <cmath>
#include <deque>

int to_index(deque<unsigned int> q) {
    unsigned int index = 0;

    for(unsigned int i = 0; i < q.size(); ++i) {
        index += (pow(3, i) * q[i]);
    }

    return index;
}

int main() {
    deque<unsigned int> q;
    unsigned int num_alleles = 3;
    unsigned int ndim = 3;
    unsigned int tmp;
    unsigned int tmp2;

    unsigned int j = 0;
    
    for(unsigned int i = 0; i < ndim; ++i) {
        q.push_front(0);
    }
    
    while(not q.empty()) {
        
        if (q.size() == ndim) {
            printf("* ");
            for(unsigned int i = 0; i < q.size(); ++i) {
                printf("%d ", q[i]);
            }
            printf(" (%d) ?= %d\n", (int) j++, to_index(q));
            
        }
            
        tmp = q.front() + 1;
        q.pop_front();
        
        if(tmp < num_alleles) {
            q.push_front(tmp);
            tmp2 = ndim - q.size();
            for(unsigned int i = 0; i < tmp2; ++i) {
                q.push_front(0);
            }
        }
    }    
}
