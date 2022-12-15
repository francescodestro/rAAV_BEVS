# Mechanistic model for rAAV manufacturing in the baculovirus expression vector system 

v1: batch model (synchronous infection, up to two waves) for TwoBac (Smith et al., 2009) and ThreeBac (Urabe et al., 2002) systems. The model simulates baculovirus binding, baculovirus transport to nucleus and replication, release of budded baculovirus, transcription and translation of AAV genes, rAAV capsid formation, Rep protein synthesis, transgene rescue and amplification, and transgene encapsidation.

Files:
- `run_simulation.m`: external wrapper to be used for setting inputs and running simulation
- `BEVS_simulation.m`: function called by external wrapper to run simulation
- `BEVS_parameters.m`: function containing model parameters
- `BEVS_model.m`: function containing model equations

Input: 
- System type (TwoBac/ThreeBac)
- Cell density at time of infection
- Multeplicity of infection
- Simulation duration

Output: 
- State vector at time intervals of 0.1 h

How to cite: Destro, F., Barone, P. W., Srinivasan, P., Springs, S. L., Cecchini, S., Kanter, J.M., S. L., Kotin, R. M., Sinskey, A. J. and Braatz, R. D. (2023). Mechanistic modeling explains the production dynamics of recombinant adeno-associated virus with the baculovirus expression vector system. Submitted.

## State vector legend
The legend of the state vector for TwoBac and ThreeBac simulations follows. The indexing of the state vector is organized in modules consistent between TwoBac and ThreeBac, to access the variables in the output array more easily. As a result, several entries of the state vector are null, especially for TwoBac. More computationally efficient implementations can be developed, with a more complex indexing of the state vector.

### TwoBac
Index j refers to: <br>
1: cells infected by repcapBV <br>
3: cells infected by goiBV <br>
5: cells infected by repcapBV and goiBV <br>

<strong> Cells concentration </strong>  <br>
<br>
`x(1)` = uninfected viable cells density (#/mL) <br>
`x(4+j)`= concentration of viable j - first wave [#/mL]  <br>
`x(28+(j-1)*22)` = concentration of nonviable j - first wave [#/mL]  
`x(204+j)` = concentration of viable j - second wave [#/mL]  <br>
`x(228+(j-1)*22)` = concentration of nonviable j - second wave [#/mL]  

<strong> Virions </strong>  <br>
<br>
`x(2)` = repcapBV virion concentration [#/mL]  <br>
`x(4)` = goiBV virion concentration [#/mL] 

<strong>  Intracellular species </strong>  <br>
Expressed as total concentration in the system [#/mL], calculated as concentration per cell multiplied by cell density <br>
Indexes are reported for the first wave. For the second wave, add 200 to each index. <br>
<br>
`x(12+(j-1)*22)` = repcapBV bound to receptors of viable j [#/mL]  <br>
`x(14+(j-1)*22)` = goiBV bound to receptors of viable j [#/mL]  <br>
`x(15+(j-1)*22)` = repcapBV DNA in nucleus of viable j [#/mL] <br>
`x(17+(j-1)*22)` = goiBV DNA in nucleus of viable j [#/mL] <br>
`x(19+(j-1)*22)` = rep mRNA in viable j [#/mL]  <br>
`x(20+(j-1)*22)` = cap mRNA in viable j [#/mL]  <br>
`x(21+(j-1)*22)` = transgene mRNA in viable j [#/mL]  <br>
`x(22+(j-1)*22)` = transgene copies in viable j [#/mL]  <br>
`x(23+(j-1)*22)` = Rep78 in viable j [#/mL]  <br>
`x(24+(j-1)*22)` = Rep52 in viable j [#/mL]  <br>
`x(25+(j-1)*22)` = empty rAAV capsids in viable j [#/mL]  <br>
`x(26+(j-1)*20)` = GOI protein in viable j [#/mL]  <br>
`x(27+(j-1)*22)` = filled rAAV capsids in viable j [#/mL]  <br>
`x(29+(j-1)*22)` = transgene copies in nonviable j [#/mL]  <br>
`x(30+(j-1)*22)` = Rep78 in nonviable j [#/mL]  <br>
`x(31+(j-1)*22)` = Rep52 in nonviable j [#/mL]  <br>
`x(32+(j-1)*22)` = empty rAAV capsids in nonviable j [#/mL]  <br>
`x(33+(j-1)*22)` = filled rAAV capsids in nonviable j [#/mL]  <br>

### ThreeBac
Index j refers to: <br>
1: cells infected by repBV <br>
2: cells infected by capBV <br>
3: cells infected by goiBV <br>
4: cells infected by repBV and capBV <br>
5: cells infected by repBV and goiBV <br>
6: cells infected by goiBV and capBV <br>
7: cells infected by repBV, capBV and goiBV  <br>

<strong> Cells concentration </strong>  <br>
<br>
`x(1)` = uninfected viable cells density (#/mL) <br>
`x(4+j)`= concentration of viable j - first wave [#/mL]  <br>
`x(28+(j-1)*22)` = concentration of nonviable j - first wave [#/mL]  
`x(204+j)` = concentration of viable j - second wave [#/mL]  <br>
`x(228+(j-1)*22)` = concentration of nonviable j - second wave [#/mL]  

<strong> Virions </strong>  <br>
<br>
`x(2)` = repBV virion concentration [#/mL]  <br>
`x(3)` = capBV virion concentration [#/mL]  <br>
`x(4)` = goiBV virion concentration [#/mL] 

<strong>  Intracellular species </strong>  <br>
Expressed as total concentration in the system [#/mL], calculated as concentration per cell multiplied by cell density <br>
Indexes are reported for the first wave. For the second wave, add 200 to each index. <br>
<br>
`x(12+(j-1)*22)` = repBV bound to receptors of viable j [#/mL]  <br>
`x(13+(j-1)*22)` = capBV bound to receptors of viable j [#/mL]  <br>
`x(14+(j-1)*22)` = goiBV bound to receptors of viable j [#/mL]  <br>
`x(15+(j-1)*22)` = repBV DNA in nucleus of viable j [#/mL] <br>
`x(16+(j-1)*22)` = capBV DNA in nucleus of viable j [#/mL] <br>
`x(17+(j-1)*22)` = goiBV DNA in nucleus of viable j [#/mL] <br>
`x(18+(j-1)*22)` = rep78 mRNA in viable j [#/mL]  <br>
`x(19+(j-1)*22)` = rep52 mRNA in viable j [#/mL]  <br>
`x(20+(j-1)*22)` = cap mRNA in viable j [#/mL]  <br>
`x(21+(j-1)*22)` = transgene mRNA in viable j [#/mL]  <br>
`x(22+(j-1)*22)` = transgene copies in viable j [#/mL]  <br>
`x(23+(j-1)*22)` = Rep78 in viable j [#/mL]  <br>
`x(24+(j-1)*22)` = Rep52 in viable j [#/mL]  <br>
`x(25+(j-1)*22)` = empty rAAV capsids in viable j [#/mL]  <br>
`x(26+(j-1)*20)` = GOI protein in viable j [#/mL]  <br>
`x(27+(j-1)*22)` = filled rAAV capsids in viable j [#/mL]  <br>
`x(29+(j-1)*22)` = transgene copies in nonviable j [#/mL]  <br>
`x(30+(j-1)*22)` = Rep78 in nonviable j [#/mL]  <br>
`x(31+(j-1)*22)` = Rep52 in nonviable j [#/mL]  <br>
`x(32+(j-1)*22)` = empty rAAV capsids in nonviable j [#/mL]  <br>
`x(33+(j-1)*22)` = filled rAAV capsids in nonviable j [#/mL]  <br>
