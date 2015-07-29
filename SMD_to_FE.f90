
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!@                                                                   @
!@ This program uses Jarzynski relationship to calculate free energy @
!@ profiles from work done in independent SMD simulations            @
!@                                                                   @
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
! Program written by Jason Swails, jason.swails@gmail.com
! Last update 7/23/2008
!
!
! Description: This program uses the Jarzynski relationship to convert a collection of
! work data compiled from multiple steered molecular dynamics simulations into a free
! energy profile for the process.  I have a script that complements this program,
! compiling a single file in the format read by this program from raw data files generated
! by AMBER's SMD module.

program SMD_to_FE
implicit none

real, dimension(:), allocatable :: Work_array

real, parameter :: KB = 1.98722E-3
real :: TEMP, distance, BETA, work, holder
integer :: num_trials
character (len=40) :: dataFile, outputFile
character (len=3) :: response
character (len=1) :: input_open
integer :: number_steps, i, c

! Work_array stores the values of work for each simulation: updated every new distance
! KB = Boltzmann's constant in kcal/mol*K (R, or gas constant)
! TEMP = system temperature, input by the user and checked for mistakes
! distance = the distance at which the Work values occur in the file
! num_trials = number of trials
! dataFile = input name for the file containing the data
! outputFile = name of file created to store output
! response = yes/no response to queries
! input_open = indicator so that error message is not displayed before error is made
! number_steps = counter that keeps track of the number of steps written to jar.log by AMBER
! i, c = counters used in loops to keep track of indices
! BETA = stat-mech beta, 1/kT
! work = keeps track of exponential of work 
! holder = holding variable that holds a variable's value while it is changed

!==================================================================!

! Get the temperature and impose checks to make sure it was entered correctly, i.e. Kelvins
! and no obvious typos made

!============================ GET TEMPERATURE =============================================
do
    print*, "Enter the temperature of the system."
    read (*,*) TEMP

    if (TEMP < 0) then
        print*, "Invalid Temperature.  Enter system temp in K"
    else if (TEMP < 100) then
        write (*,*) "Is the temperature really ", TEMP, " K?"
        read (*,*) response
        if (response == "yes" .or. response == "Yes" .or. response == "YES" .or. response == "y" &
           .or. response == "Y" .or. response == "YEs" .or. response == "yES") then
            exit
        end if ! response == ...
    else if (TEMP > 1000) then
        write (*,*) "Is the temperature really ", TEMP, " K?"
        read (*,*) response
        if (response == "yes" .or. response == "Yes" .or. response == "YES" .or. response == "y" &
           .or. response == "Y" .or. response == "YEs" .or. response == "yES") then
            exit
        end if ! response == ...
    else
        exit
    end if ! TEMP <, etc.
end do
!========================= END OF GET TEMPERATURE ========================================

! Ask for the input file, and check to make sure it exists.  User can exit by typing
! quit here or in the output file prompt, which comes right after

input_open = "n"

10 if ( input_open == "y" ) then
    print*, "There was an error opening the file, enter another one."
end if

print*,"Enter the name of the file with the data."
read(*,*) dataFile
if ( dataFile == "quit" ) then
    goto 12
end if

input_open = "y"
open(7,file=dataFile,STATUS='OLD',ERR=10)

print*,"Enter the desired name of the output file."
read(*,*) outputFile
if( outputFile == "quit" ) then
    goto 12
end if


!======================== END OF GET FILE NAMES ========================================

! Here, determine the number of trials that were done
num_trials = 0
i = 1
do
    num_trials = num_trials + 1
    read(7,FMT='(F8.5)',eor=100,advance='NO') holder
    i = i + 1
end do
100 rewind(7)

!======================== END OF TRIAL NUM DETERMINATION ===============================
! Set initial values
work = 0
BETA = 1 / (KB * TEMP)
number_steps = 1
allocate (Work_array(num_trials))

! Create the output file
open(8, file=outputFile)


! Implement Jarzynski's method
! exp(FE(x)/kT) = < exp(W(x)/kT) >
do
    read(7,*,end=11)(Work_array(i), i=1,num_trials)
    distance = Work_array(1)

    do c=2, num_trials
        work = exp(-BETA * Work_array(c)) + work
    end do

    holder = - log (work / (num_trials - 1) ) / BETA
    write (8,*) distance, holder
    number_steps = number_steps + 1
    work = 0
end do


! End the program, deallocate the array and close the files
11 continue
deallocate(Work_array)
close(7)
close(8)
print*, "Your free energy profile is complete! ", outputFile

12 continue
if( dataFile == "quit" .or. outputFile == "quit") then
    print*, "Program terminated.  Nothing computed."
end if

end program SMD_to_FE
