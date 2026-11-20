# Credits

## Authors
**Adam Ruiz Herrero** and **Ali Yamani El Yamani**, joint and equal authors of the whole
project: lab sessions 1–4 and the 2026 extensions (marked **[2026]** in the code and the
README). The exploration sections of `scripts/s04_pulse_compression.m` come from study demos
Adam wrote for Lab 4.

## Course
*Radiolocalización*, UPC-EETAC (Castelldefels). The lab statements specify the inputs and outputs
of each function, and are not redistributed here. The only code they include is two RCS routines
in Lab 5 (sphere and circular plate), adapted by the course from B. R. Mahafza, *Radar Systems
Analysis and Design Using MATLAB*, 2005. They are not included: `src/rcs` holds independent
implementations (exact Mie series, physical optics), and the course routines were used locally
only to cross-check them.

## References
- J. I. Marcum, "A statistical theory of target detection by pulsed radar", 1960 (Marcum Q function).
- M. I. Skolnik, *Introduction to Radar Systems*: detection curves and Albersheim reference values.
- Albersheim's empirical detection equation, as given in the course formula sheet. The sheet
  printed the last coefficient as 0.17; 1.7 is the value that reproduces the reference curves
  (Pd = 0.9, Pfa = 10⁻⁶ → 13.1 dB).

## Development assistance
Claude (Anthropic) acted as a development assistant for code refactoring, test scaffolding and
documentation. The engineering design, the analysis and the authorship belong to Adam and Ali.
