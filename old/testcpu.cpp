#include <iostream>
#include <cinttypes>
#include <cmath>
#include <iomanip>
#include <vector>
#include <fstream>
using namespace std;

#define printFeature(feature) if(builtin_cpu_supports(feature)) cout << feature << endl;
void sseTest()
{
    builtin_cpu_init();
    printFeature("mmx");
    printFeature("sse");
    printFeature("sse2");
    printFeature("sse3");
    printFeature("ssse3");
    printFeature("sse4.1");
    printFeature("sse4.2");
    printFeature("avx");
    printFeature("avx2");
}

int main()

{
	sseTest();
}
