/**
 * Implements an exception for solvers which are not initialized.
 */

#ifndef EIGENVALUE_PROJECT_UNINITIALIZEDSOLVER_H
#define EIGENVALUE_PROJECT_UNINITIALIZEDSOLVER_H

class UninitializedSolver: public std::invalid_argument {
private:
    std::string mVar;
    std::string mMsg;
public:
    UninitializedSolver(std::string var, std::string msg = "") :
            std::invalid_argument(msg.c_str()) {
        mVar = var;
        mMsg = msg;
    }
    const char * what () const throw () {
        return ("UninitializedSolver: The variable " + mVar
                + " has not been initialized.\n" + mMsg).c_str();
    }
};

#endif //EIGENVALUE_PROJECT_UNINITIALIZEDSOLVER_H
