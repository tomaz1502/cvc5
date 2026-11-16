/******************************************************************************
 * This file is part of the cvc5 project.
 *
 * Copyright (c) 2009-2026 by the authors listed in the file AUTHORS
 * in the top-level source directory and their institutional affiliations.
 * All rights reserved.  See the file COPYING in the top-level source
 * directory for licensing information.
 * ****************************************************************************
 *
 * Checks whether CoCoA has been patched.
 */
#include <iostream>
#include <CoCoA/TmpGPoly.H>

int main()
{
    std::cout << "CoCoA::handlersEnabled = ";
    std::cout << CoCoA::handlersEnabled << std::endl;
    // The reduction step handler must also take the power product that the
    // step cancels. This fails to compile with a CoCoA that only has the
    // older version of the patch.
    CoCoA::reductionStepHandler = [](CoCoA::ConstRefRingElem,
                                     CoCoA::ConstRefPPMonoidElem) {};
    CoCoA::reductionStepHandler = nullptr;
    return 0;
}