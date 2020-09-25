# insilico-vaccine-design-code
In silico Vaccine Design
In-silico vaccine design takes advantage of the rapidly accumulating omics data: genomics and
proteomic data is of particular interest in vaccine design. The continued development of
bioinformatics tools has further driven the in-silico vaccine design process. Machine learning
algorithms have increased the predictive ability to detect both B and T cell epitopes for creating
vaccine constructs. Vaccine design relies on the ability to distinguish patterns on antigenic
peptides that can potentially stimulate an immune response with ones that don’t. (Rötzschke et
al., 1991) suggested that the predicted patterns could be used to find peptides that would bind to
MHC, engage T cell receptors and stimulate T cell response, peptides that can elicit this response
are T-cell epitopes. Bioinformatics tools have accelerated the prediction of these motifs and at
the same time accelerated development of epitope-based vaccines. B cell epitopes have the
ability to bind to immunoglobulins or antibodies and are of two type’s linear and conformational
also known as continuous and discontinuous epitopes. A combination of B and T cell epitopes is
often used to create chimeric peptides than can be used as vaccines. Predicted T-cell epitopes are
usually at an optimal length of 9-10 mers while B-cell epitopes having varying lengths but
usually the best being 15-22 amino acids in length.

The python code(generate_epitope_combinations.py) for the generation of multi-epitope vaccine constructs,this code is essential for the generation of vaccine constructs, 
the user has latitude to select the adjuvant and linkers for the B and T cell epitopes(edit the script to appropriate adjuvant and linkers).

usage: generate_epitope_combinations.py [-h] [-t TCELL] [-b BCELL]
                                        [-r RANDOMSIZE] [-tcl TCLINKER]
                                        [-bcl BCLINKER]

optional arguments:
  -h, --help            show this help message and exit
  -r RANDOMSIZE, --randomsize RANDOMSIZE
                        number of randomized epitopes
  -tcl TCLINKER, --tclinker TCLINKER
                        T-cell epitopes linker
  -bcl BCLINKER, --bclinker BCLINKER
                        B-cell epitopes linker

required arguments:
  -t TCELL, --tcell TCELL
                        T-cell epitope file path
  -b BCELL, --bcell BCELL
                        B-cell epitope file path
                        
 provide B & T cell epitopes in csv format as follows:                       
pos	epitope
255	SEEKDTNSEEDPEAEEDPDS
295	IIPSPKPLTPEQQQERELKL
272	SSSNGSSSSNSTSSSSSSTT
251	ITKAIKKPNSGSTTSSSSNT
24	NCKCHNNNSNSSSNNDTLGG
438	VNSVSTVSPVNPVNPVNPVV
490	AVNTSNPSNPVNTVNQVVNE


Gromacs,molecular dynamics simulation of vaccine-receptor complex script(gromacs_script.sbatch) also available, this is essential to study the stability 
of the vaccine-receptor complex insilico.
