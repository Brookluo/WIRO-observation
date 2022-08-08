# Introduction to CCD

A Charge-Couple Device (CCD) is made of an array of Metal-Oxide-Semiconductor (MOS) capacitors. It is widely used in astronomy and remote sensing to take high resolution pictures. The physics behind CCD is introduced here. 

## Basic Detection Mechanism 

When a photon hits a metal, its kinetic energy will be transferred to an one electron, and "kick out" that electron from its original position. MOS capacitors contains electrons stored inside the potential well created by the circuit. 
When a photon from the source hits the capacitor, it will give an electron inside the well to jump out of the well it was in. Please see the figure 1. Then, this capacitor will carry one positive charge because the electron left created a hole. A CCD is made of the MOS capacitors. An array of CCD is called parallel registers. 

<img src="./IMG_D2E4DC37CE3A-1.jpeg" alt="figure 1" width="400"/>

Figure 1

## Basic Readout Mechanism

After detection, to reconstruct the image we just captured, we need to count how many photons have hit the detector. This step is called readout. To readout a CCD pixel, the parallel registers will be moved to a serial register section. Charges will one CCD will shift to the next one. The charge on the MOS will recouple so how mange charges that a MOS carries can be counted during the coupling process. Then an output amplifier will increase the signal the process created. Since this process creates analog signals and modern computers store information in digital format, an Analog-to-Digital Converter (ADC) or Analog-to-Digital Unit (ADU) is used to transform this signal from analog signal to digital and then store the information in computer. CCD will readout one by one, row by row. Therefore, it is called serial register. As shown in figure 2. 

<img src="./Image.jpeg" alt="figure 1" width="400"/>

Figure 2

## Noises and Operation

Typically CCDs have noises from the detection steps, readout steps and the physical properties by CCD itself. In detection, because of quantum effect, there are some photons might not be able to kick out an electron. In readout step, the ADC itself has some uncertainties so the readout is not always accurate. 
Due to properties of MOS, which has a insulation layer, but the circuit itself always has some voltage to operate. Then MOS is not perfectly insulated, so there is some current flow through that is not caused by the photon from the source. This is called dark current. CCDs need to be cooled enough to reduce these noises. 

