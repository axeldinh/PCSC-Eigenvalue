//
// Created by axeld on 01/12/2021.
//

/**
 * Implements an exception for solvers which are not initialized.
 */

#ifndef EIGENVALUE_PROJECT_UNINITIALIZEDSOLVER_H
#define EIGENVALUE_PROJECT_UNINITIALIZEDSOLVER_H

using namespace std;

class UninitializedSolver: public invalid_argument {
private:
    string mVar;
    string mMsg;
public:
    UninitializedSolver(string var, string msg = "") :
            invalid_argument(msg.c_str()) {
        mVar = var;
        mMsg = msg;
    }
    const char * what () const throw () {
        return ("UninitializedSolver: The variable " + mVar
                + " has not been initialized.\n" + mMsg).c_str();
    }
};

#endif //EIGENVALUE_PROJECT_UNINITIALIZEDSOLVER_H
