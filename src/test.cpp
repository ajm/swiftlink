using namespace std;

#include <cstdio>
#include <deque>

int main() {
    deque<unsigned int> q;
    unsigned int num_alleles = 3;
    unsigned int ndim = 4;
    unsigned int tmp;
    unsigned int tmp2;
    
    for(unsigned int i = 0; i < ndim; ++i) {
        q.push_front(0);
    }
    
    while(not q.empty()) {
        
        if (q.size() == ndim) {
            printf("* ");
            for(unsigned int i = 0; i < q.size(); ++i) {
                printf("%d ", q[i]);
            }
            printf("\n");
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
