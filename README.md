# parallel-radix-sort
Parallelization of the Radix Sort algorithm using the programming language C and parallel proccessing standards.
This code was developed as a project assignment for the undergraduate Software & Programming of High Perfomance Systems course of [CEID](wwww.ceid.upatras.gr) by Filopoimin Lykokanellos. The implementation is based on the enclosed algorithms in <i>Dani Jiménez-González, Josep-L. Larriba-Pey, Juan J. Navarro.“Case Study: Memory Conscious Parallel Sorting”, Chapter 16, Algorithms for Memory Hierarchies: Advanced Lectures, pp 355-377, vol 2625, 2003</i>.
<p>
Besides the sequential algorithm, three parallel implementations were developed. 
<ol>
<li>
The shared memory paradigm was applied to the first version through the OpenMP&reg API. 
</li>
<li>
The second implementation exploits the Open MPI library project as the message passing paradigm.  
</li>
<li>
Finally a hybrid model was used by combining the two previous standards.  
</li>
</ol>
</p>
