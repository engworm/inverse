#include <iostream>
#include <string>

enum traveller_type {gomi, gomigomi};

int main(int argc, char *argv[]) {
  char *who = argv[1];


  switch (who[0]) {
    case 'p' : std::cout << "potential" << std::endl; break;
    case 'g' : std::cout << "gravity" << std::endl; break;
    default: break; 
  }
  return 0;
}
