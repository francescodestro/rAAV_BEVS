# Mechanistic model for rAAV manufacturing with the baculovirus expression vector system 

Batch model (synchronous infection, up to two waves) for TwoBac (Smith et al., 2009) and ThreeBac (Urabe et al., 2002) systems. 

This code simulates baculovirus infection and rAAV maturation in insect cells. The following phenomena are considered in the model: baculovirus binding, transport to nucleus and replication, release of budded baculovirus, transcription and translation of AAV genes, rAAV capsid formation, Rep protein synthesis, transgene rescue and amplification, and transgene encapsidation.

Files:
- `run_simulation.m`: external wrapper to be used for setting inputs and running simulation
- `BEVS_simulation.m`: function called by external wrapper to run simulation
- `BEVS_parameters.m`: function containing model parameters
- `BEVS_model.m`: function containing model equations

Input: 
- System type (TwoBac/ThreeBac)
- Cell density at time of infection
- Multiplicity of infection
- Simulation duration

Output: 
- State vector at time intervals of 0.1 h

How to cite: Destro, F., Barone, P. W., Srinivasan, P., Springs, S. L., Cecchini, S., Kanter, J.M., S. L., Kotin, R. M., Sinskey, A. J. and Braatz, R. D. (2023). Mechanistic modeling explains the production dynamics of recombinant adeno-associated virus with the baculovirus expression vector system. _Mol Ther Methods Clin Dev_ In press.

## License
The code in this repository is provided under a CC BY-NC-ND 4.0 license, as detailed in the `LICENSE` file.

## State vector legend
This section outlines the legend of the state vector for TwoBac and ThreeBac simulations. The indexing of the state vector is organized in modules consistent between TwoBac and ThreeBac, to access the variables in the output array more easily. As a result, several entries of the state vector are null, especially for TwoBac.

### TwoBac
Index `i` refers to: <br>
`i=1`: cells infected by repcapBV <br>
`i=3`: cells infected by goiBV <br>
`i=5`: cells infected by repcapBV and goiBV <br>

<strong> Cells concentration </strong>  <br>
<br>
`x(1)` = uninfected viable cells density (#/mL) <br>
`x(4+i)`= concentration of viable i - first wave [#/mL]  <br>
`x(28+(i-1)*22)` = concentration of nonviable i - first wave [#/mL]  
`x(204+i)` = concentration of viable i - second wave [#/mL]  <br>
`x(228+(i-1)*22)` = concentration of nonviable i - second wave [#/mL]  

<strong> Virions </strong>  <br>
<br>
`x(2)` = repcapBV virion concentration [#/mL]  <br>
`x(4)` = goiBV virion concentration [#/mL] 

<strong>  Intracellular species </strong>  <br>
Expressed as total concentration in the system [#/mL], calculated as concentration per cell multiplied by cell density <br>
Indexes are reported for the first wave. For the second wave, add 200 to each index. <br>
<br>
`x(12+(i-1)*22)` = repcapBV bound to receptors of viable i [#/mL]  <br>
`x(14+(i-1)*22)` = goiBV bound to receptors of viable i [#/mL]  <br>
`x(15+(i-1)*22)` = repcapBV DNA in nucleus of viable i [#/mL] <br>
`x(17+(i-1)*22)` = goiBV DNA in nucleus of viable i [#/mL] <br>
`x(19+(i-1)*22)` = rep mRNA in viable i [#/mL]  <br>
`x(20+(i-1)*22)` = cap mRNA in viable i [#/mL]  <br>
`x(21+(i-1)*22)` = transgene mRNA in viable i [#/mL]  <br>
`x(22+(i-1)*22)` = transgene copies in viable i [#/mL]  <br>
`x(23+(i-1)*22)` = Rep78 in viable i [#/mL]  <br>
`x(24+(i-1)*22)` = Rep52 in viable i [#/mL]  <br>
`x(25+(i-1)*22)` = empty rAAV capsids in viable i [#/mL]  <br>
`x(26+(i-1)*20)` = GOI protein in viable i [#/mL]  <br>
`x(27+(i-1)*22)` = filled rAAV capsids in viable i [#/mL]  <br>
`x(29+(i-1)*22)` = transgene copies in nonviable i [#/mL]  <br>
`x(30+(i-1)*22)` = Rep78 in nonviable i [#/mL]  <br>
`x(31+(i-1)*22)` = Rep52 in nonviable i [#/mL]  <br>
`x(32+(i-1)*22)` = empty rAAV capsids in nonviable i [#/mL]  <br>
`x(33+(i-1)*22)` = filled rAAV capsids in nonviable i [#/mL]  <br>

### ThreeBac
Index `i` refers to: <br>
`i=1`: cells infected by repBV <br>
`i=2`: cells infected by capBV <br>
`i=3`: cells infected by goiBV <br>
`i=4`: cells infected by repBV and capBV <br>
`i=5`: cells infected by repBV and goiBV <br>
`i=6`: cells infected by goiBV and capBV <br>
`i=7`: cells infected by repBV, capBV and goiBV  <br>

<strong> Cells concentration </strong>  <br>
<br>
`x(1)` = uninfected viable cells density (#/mL) <br>
`x(4+i)`= concentration of viable i - first wave [#/mL]  <br>
`x(28+(i-1)*22)` = concentration of nonviable i - first wave [#/mL]  
`x(204+i)` = concentration of viable i - second wave [#/mL]  <br>
`x(228+(i-1)*22)` = concentration of nonviable i - second wave [#/mL]  

<strong> Virions </strong>  <br>
<br>
`x(2)` = repBV virion concentration [#/mL]  <br>
`x(3)` = capBV virion concentration [#/mL]  <br>
`x(4)` = goiBV virion concentration [#/mL] 

<strong>  Intracellular species </strong>  <br>
Expressed as total concentration in the system [#/mL], calculated as concentration per cell multiplied by cell density <br>
Indexes are reported for the first wave. For the second wave, add 200 to each index. <br>
<br>
`x(12+(i-1)*22)` = repBV bound to receptors of viable i [#/mL]  <br>
`x(13+(i-1)*22)` = capBV bound to receptors of viable i [#/mL]  <br>
`x(14+(i-1)*22)` = goiBV bound to receptors of viable i [#/mL]  <br>
`x(15+(i-1)*22)` = repBV DNA in nucleus of viable i [#/mL] <br>
`x(16+(i-1)*22)` = capBV DNA in nucleus of viable i [#/mL] <br>
`x(17+(i-1)*22)` = goiBV DNA in nucleus of viable i [#/mL] <br>
`x(18+(i-1)*22)` = rep78 mRNA in viable i [#/mL]  <br>
`x(19+(i-1)*22)` = rep52 mRNA in viable i [#/mL]  <br>
`x(20+(i-1)*22)` = cap mRNA in viable i [#/mL]  <br>
`x(21+(i-1)*22)` = transgene mRNA in viable i [#/mL]  <br>
`x(22+(i-1)*22)` = transgene copies in viable i [#/mL]  <br>
`x(23+(i-1)*22)` = Rep78 in viable i [#/mL]  <br>
`x(24+(i-1)*22)` = Rep52 in viable i [#/mL]  <br>
`x(25+(i-1)*22)` = empty rAAV capsids in viable i [#/mL]  <br>
`x(26+(i-1)*20)` = GOI protein in viable i [#/mL]  <br>
`x(27+(i-1)*22)` = filled rAAV capsids in viable i [#/mL]  <br>
`x(29+(i-1)*22)` = transgene copies in nonviable i [#/mL]  <br>
`x(30+(i-1)*22)` = Rep78 in nonviable i [#/mL]  <br>
`x(31+(i-1)*22)` = Rep52 in nonviable i [#/mL]  <br>
`x(32+(i-1)*22)` = empty rAAV capsids in nonviable i [#/mL]  <br>
`x(33+(i-1)*22)` = filled rAAV capsids in nonviable i [#/mL]  <br>
