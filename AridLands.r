# Spatially explitic model for arid vegetation on a flat surface
# The model follows the paper by Rietkerk et al, AmNat 160, October 2002.
#
# Rietkerk, M., Boerlijst, M.C., Van Langevelde, F., HilleRisLambers, R., 
# Van de Koppel, J., Kumar, L., Klausmeier, C.A., Prins, H.H.T., De Roos, A. 2002.
# The origin of pattern formation in arid vegetation. The American Naturalist 160:524-530.
# 
# The model has been changed slightly (different diffusion coefficients) to allow  
# a fast simulation.

remove(list=ls())     # Remove all variables from memory
on=1;off=0 

require("fields")

setwd('/Users/johank/Simulations/AridLands')

Slope=off

# Parameters            Original value    Explanation and Units
R       =   0.8         # 1.00              rainfall(mm.d-1), original values 0.75, 1.00, and 1.25 
alpha   =   0.2         # 0.2               prop of surface water (d-1)
W0      =   0.2         # 0.2               bare soil infiltration rate (-)
rw      =   0.2         # 0.2               Soil water loss rate due to seepage and evaporation (d-1)
c       =    10         # 10                Plant uptake constant (g.mm-1.m-2)
gmax    =  0.05         # 0.05              Plant growth constant (mm.g-1.m2.d-1)   
d       =  0.25         # 0.25              Plant senescence rate (d-1)
k1      =     5         # 5                 Half saturation constant for plant uptake and growth (mm)
k2      =     5         # 5                 Half saturation constant for water infiltration (g.m-2)

DifP=0.1                # 0.1               The diffusion constant of P
DifW=0.1                # 0.1               The diffusion constant of W
DifO=100                # 100               The diffusion constant of O
AdvO=10                 # 10                The advection constant for water running downslopw

# Simulation parameters
Length = 32            # Size of the simulated landscape
m=64                    # Number of gridcells along the axis
NX=m                    # Number of gridcells in the x-direction
NY=m                    # Number of gridcells in the y-direction

DeltaX=Length/NX        # 8 The x dimension of a cell, in meters
DeltaY=Length/NX        # 8 The y dimension of a cell, in meters

dT=0.0005                # Timestep
Time=1                  # Begin time       
EndTime=500             # End time
NoFrames=500            # The number of times the plot is updated

frac=0.10               # The fraction of the cells that is initially filled with plants
mp  =  m+1              # mp is used in specific array designations

WinWidth=12             #      - Width of the simulation window 
WinHeight=4             #      - Height of the simulation window

# ------ Function definitions
Laplacian <- function (C, dx, dy) {
  
  NX=dim(C)[2]
  NY=dim(C)[1]
  
  dC=matrix(nrow=NY, ncol=NX, data=0)
  
  dC[2:(NY-1),2:(NX-1)] = 
    (C[2:(NY-1),1:(NX-2)] + C[2:(NY-1),3:NX] - 2*C[2:(NY-1),2:(NX-1)])/dx/dx +
    (C[1:(NY-2),2:(NX-1)] + C[3:NY,2:(NX-1)] - 2*C[2:(NY-1),2:(NX-1)])/dy/dy
  
  return(dC)
  
}

GradientX <- function (C, dx) {
  
  NX=dim(C)[2]
  NY=dim(C)[1]
  
  dC=matrix(nrow=NY, ncol=NX, data=0)
  
  dC[2:(NY-1),2:(NX-1)] = 
    - (C[2:(NY-1),2:(NX-1)] - C[2:(NY-1),3:NX] )/dx/2
  
  return(dC)
  
}

PeriodicBoundary <- function (C) {
  
  NX=dim(C)[2]
  NY=dim(C)[1]
  
  C_new=C
  
  C_new[1,]=C_new[m-1,] 
  C_new[m,]=C_new[2,] 
  C_new[,1]=C_new[,m-1] 
  C_new[,m]=C_new[,2]
  
  return(C_new)
  
}

# Initialisation: declaring variable and setting up initial values
# All arrays of dimension m x m
P = W = O = dP = dW = dO = NetP = NetW = NetO = matrix(nrow=m,ncol=m) 

# A an array to store time and vegetation biomass
TimeRec=VegRec=vector(length=NoFrames) 

#------ Initial setup and calculation ----------------------------------

# initial state of A and M at the start of the run
P=+matrix(ncol=m,nrow=m,data=(runif(m*m)<frac)*25)  # At R=1.2, for Symmetry breaking
W[,]=R/rw          # Equilibrium water level at full water infiltration divided by 2 
O[,]=R/alpha*2     # Equilibrium water level at zero plant biomass 

DeltaX=Length/m    # The x dimension of a cell, in meters
DeltaY=Length/m    # The y dimension of a cell, in meters

Time =  0          # Begin time 
ii   =  1e6        # Setting the plot counter to max, so that drawing start immediately
jj   =  0          # The counter needed for recording data during the run

# ------- Setting up the figure ------------------------------------------
## Open a graphics window (Darwin stands for a Mac computer)
if (Sys.info()["sysname"]=="Darwin"){
  quartz(width=WinWidth, height=WinHeight, 
         title="Arid bushland model")
} else
  windows(width = WinWidth, height = WinHeight,
          title="Arid bushland model")

par(mfrow=c(1,3), mar=c(3, 2, 2, 5) + 0.1)
  
OGraphMax = R/(alpha*W0)   # The maximal value that the O plot accepts
WGraphMax = R/rw*1.5       # The maximal value that the W plot accepts
PGraphMax = 10             # The maximal value that the P plot accepts

Bush.colors= function(x)rev(terrain.colors(x)) # The white to green color range
water.palette = colorRampPalette(c("white", "blue"))

# ------------ The simulation loop ---------------------------------------
print(system.time(       # This function calculates how long the simulation takes
while (Time<=EndTime){   # Here the time loop starts   
  
  # The rates of change due to growth/mortality/uptake/infiltration etc are calculated
  # Diffusive processes are not yet included
    
  # dO/dt = Rainfall - Infiltration
    drO   = R        - alpha*(P+k2*W0)/(P+k2)*O;
      
  # dW/dt = Infiltration                      - Uptake                   - Drainage
    drW   = alpha*(P+k2*W0)/(P+k2)*O - gmax*W/(W+k1)*P - rw*W;
      
  # dP/dt = Plant growth               - Mortality
    drP   = c*gmax*W/(W+k1)*P - d*P
  
  # Calculate Advective and Diffusive flow of water and plants

  if (Slope==on) NetO = AdvO*GradientX(O,DeltaX)
      else NetO=DifO*Laplacian(O,DeltaX,DeltaY)
  NetW=DifW*Laplacian(W,DeltaX,DeltaY)
  NetP=DifP*Laplacian(P,DeltaX,DeltaY)
      
  # Summing up local processes and lateral flow to calculate new A and M
  O=O+(NetO+drO)*dT 
  W=W+(NetW+drW)*dT 
  P=P+(NetP+drP)*dT  
    
  # Periodic Boundary conditions
  P = PeriodicBoundary(P)  
  W = PeriodicBoundary(W)  
  O = PeriodicBoundary(O)
    
  # Graphic representation of the model every now and then
  if (ii>=EndTime/NoFrames/dT)
      {# First plot, containing surface water
       image.plot(pmin(O,OGraphMax), zlim=c(0,OGraphMax), xaxt="n", yaxt="n",
             col = water.palette(255),asp=1, bty="n",
             legend.shrink = 0.95, legend.width = 1.8)
       title("Surface water", line=0.5)      

       # Second plot, containing soil water
       image.plot(pmin(W,WGraphMax), zlim=c(0,WGraphMax), xaxt="n", yaxt="n",
             col = water.palette(255),asp=1, bty="n",
             legend.shrink = 0.95, legend.width = 1.8)       
       title("Soil water", 
              xlab = paste("Time : ",sprintf("%1.0f",Time),
                     "of" ,sprintf("%1.0f",EndTime), "days"), 
                     cex=1.5, line=0.5)

       # Last plot, containing vegetation biomass
       image.plot(pmin(P,PGraphMax), zlim=c(0,PGraphMax), xaxt="n", yaxt="n",
             col = Bush.colors(255),asp=1, bty="n",
             legend.shrink = 0.95, legend.width = 1.8)
       title("Vegetation biomass", line=0.5)
       
       # The following two lines prevent flickering of the screen
       dev.flush() # Force the model to update the graphs
       dev.hold()  # Put all updating on hold
       
       ii=0    # Resetting the plot counter 
       jj=jj+1 # Increasing the Recorder counter
       
       TimeRec[jj]=Time # The time in days
       VegRec[jj]=mean(P)  # Mean mussel biomass 
      } 

  Time=Time+dT  # Incrementing time with one
  ii=ii+1       # Incrementing the plot counter with one
 
}  # Here the time loop ends

))

## Open a graphics window (Darwin stands for a Mac computer)
if (Sys.info()["sysname"]=="Darwin"){
  quartz(width=WinWidth, height=WinHeight, 
         title="Arid bushland model")
} else
  windows(width = WinWidth, height = WinHeight,
          title="Arid bushland model")

plot(TimeRec,VegRec, xlab="Time (days)", ylab="Vegetation biomass")

