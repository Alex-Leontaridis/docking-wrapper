Architecture
============

.. include:: ../../docs/PATH_AND_LOGGING_REFACTOR_SUMMARY.md

.. note::
   The following Mermaid diagram illustrates the high-level architecture:

.. mermaid::

   graph TD;
     A[User CLI/API] --> B[Batch Pipeline]
     B --> C1[Structure Preparation]
     B --> C2[Docking Engines]
     B --> C3[ML Models]
     B --> C4[Analysis Tools]
     C2 --> D1[Vina]
     C2 --> D2[GNINA]
     C2 --> D3[DiffDock]
     C3 --> D4[EquiBind]
     C3 --> D5[UMol]
     C3 --> D6[NeuralPLexer]
     C4 --> D7[Boltz2]
     C4 --> D8[PLIP]
     C4 --> D9[FPocket]
     B --> E[Results Parsing]
     E --> F[Summary Reports] 