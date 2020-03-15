with open("index.h", 'w') as f:
    f.write("""
#include <iostream>
using namespace std;
bool boolator_a() {
    cout << "A" << endl;
    return true;
}""")
