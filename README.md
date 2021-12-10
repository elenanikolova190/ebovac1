<div id="top"></div>

<!-- PROJECT NAME -->
<br />
  <h3 align="center">EBOVAC 1 Application</h3>

  <p align="center">
    An application to simulate vaccine trials for Ebola 
    <br />
    <br />
    <a href="https://github.com/othneildrew/Best-README-Template">View Demo</a>
  </p>
</div>


<!-- TABLE OF CONTENTS -->
<details>
  <summary>Table of Contents</summary>
  <ol>
    <li>
      <a href="#about-the-project">About The Project</a>
      <ul>
        <li><a href="#built-with">Built With</a></li>
      </ul>
    </li>
    <li>
      <a href="#models">Trial Models Description</a>
      <ul>
        <li><a href="#1dlp">One Dose Leaky Protection Model</a></li>
        <li><a href="#1daon">One Dose All-or-Nothing Protection Model</a></li>
        <li><a href="#2dlp">Two Dose Leaky Protection Model</a></li>
        <li><a href="#2daon">Two Dose All-or-Nothing Protection Model</a></li>
      </ul>
    </li>
    <li><a href="#tool">Tool</a></li>
   
  </ol>
</details>

<!-- ABOUT THE PROJECT -->
## About The Project

EBOVAC 1 project was started in 2016.
This work was completed as part of Work Package 4 of the project. The purpose of this work is to develop and implement a software tool for 
simulating an Ebola epidemic model that would allow the researchers and decision makers to simulate Ebola epidemics in order to plan vaccine trials. 
It should allow to determine how trial feasibility changes with time and what trial type would be potentially suitable.

Four trial models have been developed: 
* One dose vaccine leaky protection model
* One dose vaccine all-or-nothing protection model
* Two dose vaccine leaky protection model
* Two dose vaccine all-or-nothing protection model

Details for each model have been provided below.  

<p align="right">(<a href="#top">back to top</a>)</p>


### Basic Assumptions 

Regarding the trial design, the following assumptions have been made.

 * Homogenious mixing among the population
 * Randomised control trial
 
 ###  Current Status of the Tool
 
Currently, only the one dose leaky protection model results have been shown. The list of issues presents how the tool needs to be developed further and improved. 
 
 <div id="ebovac_tool", align="center">
  <img src="https://github.com/elenanikolova190/ebovac1/blob/main/images/ebovac_interface.png"
     title="Ebovac tool"
     width = "400", 
     height: auto />
</div>
 
The application allows the users to select a number of input parameters, starting with the type of vaccine protection model.

#### POPULATION

* Affected population - the total population affected by the epidemic (default 100 000 people)
* Initial number of cases - number of infections at the beginning of the simulation (default 10 cases)

#### VACCINE

* Vaccine Efficacy - projected (default 0.75) 
* Onset of Immunity protection (in days) - onset of immunity after administering the vaccine (default 14 days)
* Onset of Immunity after Second Dose (in days) - (default 10 days)

#### TRIAL

* Starting Day of trial - How many days after the beginning of the simulation would trial begin
* Number of vaccine doses - how many doses of vaccine are to be administered during the trial (for 2 dose vaccines, how many of each dose)
* Trial Duration (in days) - how long each dose would take to administer (default 30 days)
* Time gap between 1st and 2nd dose - it is assume that the first dose administration would be complete before the second round begins. Thus the time gap should be equal or bigger than the trial duration.

#### EPIDEMIC CHARACTERISTICS

* Reproductive number (R0)- (default 1.5)
* Average latent period (in days) - (default 10 days)
* Average infectious peiod (in days)" - (default 12 days)


<!-- MODELS -->
## MODELS

### One Dose Leaky Protection Model

The compartments of this model are presented in the diagram below. After vaccination, the participants are moved to Sv (Suseptible Vaccinated) and C1 (Control arm, Compartment 1) compartments.
As it is possible that some participants might be exposed to the virus prior to vaccination, these are acounted by the Evac, Ivac and Rvac compartments.
If some participants get infected during the washout period, these would be moved to the Ev -> Iv -> Rv compartments.
After the washout period, most of the participants in the trial would move to the Vp (Vaccinated protected) compartment. As time progresses, as this is a leaky model, 
some of them might get exposed and move to Ep -> Ip -> Rp compartments for the vaccinated arm of the trial. 

In the control group, as they do not receive the real vaccine, all the individuals are susceptible to infection. To be able to compare the number of infected participants 
from the control arm with those vaccinated at each stage, two compartments, C1 and  C2 are being implemented with their corresponding successive E, I and R compartments.
 More details and the differential equations for the model can be found in the specification document here.  

<div id="1dlp_diagram", align="center">
  <img src="https://github.com/elenanikolova190/ebovac1/blob/main/images/OneD_LeakyP.png"
     title="One Dose Leaky Protection Model Diagram"
     width = "400", 
     height: auto />
</div>

### One Dose All or Nothing Protection Model

For this model, the vaccinated participants are being moved to two compartments:
* those, who have developed full immunity go into compartment  Vp 
* those, for whom the vaccine has not worked are moved into compartment  Vu. They are susceptible to infection.


<div id="1daon_diagram", align="center">
  <img src="https://github.com/elenanikolova190/ebovac1/blob/main/images/OneD_A-o-NP.png"
     title="One Dose All-or-Nothing Protection Model Diagram"
     width = "400", 
     height: auto />
</div>

### Two Dose Leaky Protection Model

In the two-dose leaky protection model there are 39 compartments. The trial is divided into four stages: a) administering the first dose,
b) onset of partial immunity, c) administering of the boost, and d) onset of full protection. After administering the first dose, 
the participants move from compartment Sv to Pv at a rate k1. Later after receiving the boost, compartment B, they move to Pb when full immunity
is achieved after a period of k2. The participants in the trial who have been given the mock vaccine (the control group) follow the same paths moving from C1 to C2 
and after receiving the boost to Cb2. That is done to keep account of the number of infections that occur at each stage of the trial in both arms. 

<div id="2dlp_diagram", align="center">
  <img src="https://github.com/elenanikolova190/ebovac1/blob/main/images/TwoD_LeakyP.png"
     title="Two Dose Leaky Protection Model Diagram"
     width = "400", 
     height: auto />
</div>


### Two Dose All or Nothing Protection Model

In this case scenario there are 42 compartments altogether. 
Here, after stage 1, the participants for whom the vaccine has been effective are moved into P1, Protected, 
compartment, while the rest are in the Vu, Unprotected, compartment. It is assumed for version one of the software that 
the efficacy from the first dose 100%, while the participants in the unprotected compartment have not developed any immunity at all. 
After the boost is administered, these are moved to compartments Bp and Bu and after another period, k2 into Paband Uab respectively. 


<div id="2daonp_diagram", align="center">
  <img src="https://github.com/elenanikolova190/ebovac1/blob/main/images/TwoD_All-o-N_P.png"
     title="Two Dose All-or-Nothing Protection Model Diagram"
     width = "400", 
     height: auto />
</div>



<p align="right">(<a href="#top">back to top</a>)</p>

