#ifndef HAVE_INTEGRATION_H
#define HAVE_INTEGRATION_H

double adapt24 ( double func(double), double leftEndPt, double rightEndPt, double absAcc, double relAcc, double second_funcVal, double third_funcVal, int numOfRecursions, double* integrationError );
double adapt ( double func(double), double leftEndPt, double rightEndPt, double absAcc, double relAcc, double* integrationError );
double openQuad(double func(double), double leftEndPt, double rightEndPt, double absAcc, double relAcc, double* integrationError );
double integrate(double func(double), double leftEndPt, double rightEndPt, double absAcc, double relAcc, double* integrationError);

#endif
