#include "chroma.h"

using namespace Chroma;

int main(int argc, char *argv[]) 
{ 
  Chroma::initialize(&argc, &argv);

  CircularBuffer<int> c(3);


  for(int i=0; i < 10; i++) { 
    c.push(i);

    QDPIO::cout << "Getting at elements: " << std::endl;
    QDPIO::cout << "Size is " << c.size() << std::endl;
    QDPIO::cout << "contents: " ;
    for(int j=0; j < c.size(); j++) { 
      int ith_elem;
      c.get(j, ith_elem);
      QDPIO::cout << " " << ith_elem;
    }
    QDPIO::cout << std::endl;
  }


  // Nuke buffer
  c.reset();

  try { 
    int foo;
    c.get(0,foo);
  }
  catch (const CircularBuffer<int>::OutOfBoundsException& e) {
    QDPIO::cout << "Caught deliberate exception: " << e.error_string << std::endl;
    QDPIO::cout << "Index requested : " << e.i << std::endl;
    QDPIO::cout << "Buffer size :" << e.size << std::endl;
  }


  for(int i=11; i < 20; i++) { 
    c.push(i);

    QDPIO::cout << "Getting at elements: " << std::endl;
    QDPIO::cout << "Size is " << c.size() << std::endl;
    QDPIO::cout << "contents: " ;
    for(int j=0; j < c.size(); j++) { 
      int ith_elem;
      c.get(j, ith_elem);
      QDPIO::cout << " " << ith_elem;
    }
    QDPIO::cout << std::endl;
  }

  try { 
    int foo;
     c.get(7, foo);
  }
  catch (const CircularBuffer<int>::OutOfBoundsException& e) {
    QDPIO::cout << "Caught deliberate exception: " << e.error_string << std::endl;
    QDPIO::cout << "Index requested : " << e.i << std::endl;
    QDPIO::cout << "Buffer size :" << e.size << std::endl;
  }


  Chroma::finalize();
  exit(0);
}
