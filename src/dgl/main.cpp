#include "data.h"
#include "dglalg.h"

using namespace riot;


int main()
{
  Data D; //Daten werden gesetzt.
  Files file(5,D.Title_); //Files für Ausgabe.

  DGLSolver::Algorithm().prepareAlgorithm(&D);
  DGLSolver::Algorithm().run(&file);

  return 0;
}
