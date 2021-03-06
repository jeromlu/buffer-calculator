    module read_parameters
    ! reads and/or sets all required variables

    ! parameters
    implicit none
    double precision :: input_molarities(11,2)
    double precision :: input_pH(2), input_V(2)
    integer :: input_phos_type(2)
    double precision :: MW(20), dens(5)
    double precision :: pKa_acet, pKa_phos(3), pKa_mala(2), pKa_succ(2), pKa_tris, pKa_citr(3), pKa_argi(3)
    double precision :: pKa_hist(3), pKa_cyst(3), pKa_EDTA(4)
    double precision :: z_acet, z_phos(3), z_mala(2), z_succ(2), z_tris, z_citr(3), z_argi(3)
    double precision :: z_hist(3), z_cyst(3), z_EDTA(4)
    double precision :: A_DH, B_DH(10)
    double precision :: Diff(23)
    double precision :: ion_rad(23), DPKDT(26), pKw, a_dh_relationship(2), visc(2),KB6pi


    ! subroutines
    contains

    subroutine read_all
    ! parameters
    implicit none
    integer :: loop

    ! reading inputs
    ! opening the file
    open(unit=20, file='param_all.dat', status='old', action='read')
    do loop = 1,11
        read(20,*) input_molarities(loop,:) ! in order, we have: acetate, phosphate, succinate, tris, citrate, malate, arginine, histidine, cysteine, EDTA, NaCl
        input_molarities(loop,:) = input_molarities(loop,:)/1000d0 ! transforming into mol/L
    end do
    read(20,*) input_pH
    read(20,*) input_phos_type
    read(20,*) input_V
    close (20)

    ! set physical inputs
    ! molecular weights
    MW(1) = 60.5d0 ! AcOH
    MW(2) = 82d0 ! NaAcO- anhydre
    MW(3) = 98d0 ! H3PO4 85%
    MW(4) = 137.98d0 ! NaH2PO4 monohydrate
    MW(5) = 141.9588d0 ! Na2HPO4
    MW(6) = 163.94d0 ! Na3PO4
    MW(7) = 118.09d0 ! succinic acid
    MW(8) = 121.13d0 ! tris-base
    MW(9) = 157.59d0 ! tris HCl
    MW(10) = 210.66d0 ! arginine HCl
    MW(11) = 155.15d0 ! histidine
    MW(12) = 209.63d0 ! histidine HCl H2O
    MW(13) = 39.99d0 ! NaOH
    MW(14) = 58.44d0 ! NaCl
    MW(15) = 18.01d0 ! H2O
    MW(16) = 36.46d0 ! HCl
    MW(17) = 210.14d0 ! citric acid monohydrate
    MW(18) = 134.09d0 ! malic acid
    MW(19) = 175.63d0 ! cysteine HCL H2O
    MW(20) = 372.24d0 ! EDTA 2xH2O 2xNa

    ! densities
    dens(1) = 997.0479d0 ! H2O
    dens(2) = 1049d0 ! AcOH
    dens(3) = 1329d0 ! 10 M NaOH
    dens(4) = 1123.55d0 ! 25% HCl
    dens(5) = 1685d0 ! 85% H3PO4

    ! pKa values
    pKa_acet = 4.76d0

    pKa_phos(1) = 2.15d0
    pKa_phos(2) = 7.21d0
    pKa_phos(3) = 12.33d0

    pKa_mala(1) = 3.4d0
    pKa_mala(2) = 5.13d0

    pKa_succ(1) = 4.21d0
    pKa_succ(2) = 5.64d0

    pKa_tris = 8.07d0

    pKa_citr(1) = 3.13d0
    pKa_citr(2) = 4.76d0
    pKa_citr(3) = 6.4d0

    pKa_argi(1) = 2.09d0
    pKa_argi(2) = 9d0
    pKa_argi(3) = 12.1d0

    pKa_hist(1) = 1.7d0
    pKa_hist(2) = 9.09d0
    pKa_hist(3) = 6.04d0

    pKa_cyst(1) = 1.91d0
    pKa_cyst(2) = 10.28d0
    pKa_cyst(3) = 8.14d0

    pKa_EDTA(1) = 2d0
    pKa_EDTA(2) = 2.7d0
    pKa_EDTA(3) = 6.2d0
    pKa_EDTA(4) = 10.31d0
    
    pKw = 14d0

    ! conjugate acid charges
    z_acet = 0d0

    z_phos(1) = 0d0
    z_phos(2) = -1d0
    z_phos(3) = -2d0

    z_mala(1) = 0d0
    z_mala(2) = -1d0

    z_succ(1) = 0d0
    z_succ(2) = -1d0

    z_tris = 1d0

    z_citr(1) = 0d0
    z_citr(2) = -1d0
    z_citr(3) = -2d0

    z_argi(1) = 0d0
    z_argi(2) = 1d0
    z_argi(3) = 1d0

    z_hist(1) = 0d0
    z_hist(2) = 1d0
    z_hist(3) = 1d0

    z_cyst(1) = 0d0
    z_cyst(2) = 1d0
    z_cyst(3) = 0d0

    z_EDTA(1) = 0d0
    z_EDTA(2) = -1d0
    z_EDTA(3) = -2d0
    z_EDTA(4) = -3d0

    ! Debye Huckel parameters
    A_DH = 0.5114d0
    ! source: Values of the Constants in the Debye—Hückel Equation for Activity Coefficients
    A_DH_relationship(1) = 0.000942521d0 ! temperature dependence slope
    A_DH_relationship(2) = 0.228544871d0  ! temperature dependence intercept

    
    do loop = 1,10
        B_DH(loop) = 0.1d0 ! acetate, phosphate, succinate, tris, citrate, malate, arginine, histidine, cysteine, EDTA
    end do
    B_DH(1) = 0.09d0 ! acetate
    B_DH(2) = 0.07d0 ! phosphate

    ! and the diffusion constants! https://is.muni.cz/el/1431/podzim2016/C4020/um/pom/Ionic_Conductivity_and_Diffusion_at_Infinite_Dilution.pdf
    ! http://www.aqion.de/site/194
    Diff(1) = 1.089d-09 ! cetate
    Diff(2) = 8.46d-10 ! phosphate
    Diff(3) = 6.9d-10
    Diff(4) = 6.12d-10
    Diff(5) = 7.83d-10 ! succinate
    Diff(6) = 7.83d-10
    Diff(7) = 6.84d-10 ! tris
    Diff(8) = 6.23d-10 ! citric
    Diff(9) = 6.23d-10
    Diff(10) = 6.23d-10
    Diff(11) = 7.83d-10 ! malate
    Diff(12) = 7.83d-10
    Diff(13) = 5.45d-10 ! arg
    Diff(14) = 5.62d-10 ! his
    Diff(15) = 6.87d-10 ! cys
    Diff(16) = 5.83d-10! EDTA
    Diff(17) = 5.83d-10
    Diff(18) = 5.83d-10
    Diff(19) = 5.83d-10
    Diff(20) = 1.33d-09! Na
    Diff(21) = 2.03d-09 ! Cl
    Diff(22) = 9.31d-09 ! H
    Diff(23) = 5.27d-09 ! OH
    
    ! radius of the consi=dered ion
    ion_rad(1)=2.25206222648314d-10
    ion_rad(2)=2.89893116387724d-10
    ion_rad(3)=3.55434168788427d-10
    ion_rad(4)=4.00734602065383d-10
    ion_rad(5)=3.13217849890184d-10
    ion_rad(6)=3.13217849890184d-10
    ion_rad(7)=3.5855201237429d-10
    ion_rad(8)=3.93659031242399d-10
    ion_rad(9)=3.93659031242399d-10
    ion_rad(10)=3.93659031242399d-10
    ion_rad(11)=3.13217849890184d-10
    ion_rad(12)=3.13217849890184d-10
    ion_rad(13)=4.49999222869751d-10
    ion_rad(14)=4.36387146733122d-10
    ion_rad(15)=3.56986283062612d-10
    ion_rad(16)=4.20668227211002d-10
    ion_rad(17)=4.20668227211002d-10
    ion_rad(18)=4.20668227211002d-10
    ion_rad(19)=4.20668227211002d-10
    ion_rad(20)=1.84398177792492d-10
    ion_rad(21)=1.20812599243357d-10
    ion_rad(22)=2.63425968274988d-11
    ion_rad(23)=4.65369215301735d-11

    ! kB/6pi
    kB6pi = 7.32d-25

    ! viscosity of water source: https://wiki.anton-paar.com/en/water/
    ! only valid between 15&35 C
    visc(1) = -2.06995d-05 ! slope 
    visc(2) = 0.007072374 ! intercept

    
    ! pKa dependency on temperature 
    ! (from: (many sources)
    ! Basic principles of electrolyte chemistry for microfluidic electrokinetics. Part I: Acid–base equilibria and pH buffers
    ! http://www.reachdevices.com/Protein/BiologicalBuffers.html
    !Temperature Dependence of the Dissociation Constants of Several Amino Acids Hidetada
    ! Computer Simulation of the Effect of Temperature on pH JAMES
    ! Biorelevant pKa (37°C) predicted from the 2D structure of the molecule and its pKa at 25°C
    dpKdT(1) = -0.0002d0! acetate
    dpKdT(2) = 0.0044d0! phosphate
    dpKdT(3) = -0.0028d0
    dpKdT(4) = -0.026d0
    dpKdT(5) = 0.000833333d0! malate
    dpKdT(6) = 0.0075d0
    dpKdT(7) = -0.0018d0 !succinate
    dpKdT(8) = 0.0002d0
    dpKdT(9) = -0.028d0! tris
    dpKdT(10) = -0.0024d0 ! citrate
    dpKdT(11) = -0.001388373d0
    dpKdT(12) = 0.001914997d0
    dpKdT(13) = 0d0! argi
    dpKdT(14) = -0.021214612d0
    dpKdT(15) = -0.005141553d0
    dpKdT(16) = 0.001296804d0! hist
    dpKdT(17) = -0.020721461d0
    dpKdT(18) = -0.012630137d0
    dpKdT(19) = 0d0 ! cyste not available
    dpKdT(20) = 0d0
    dpKdT(21) = 0d0
    dpKdT(22) = 0.009166667d0! EDTA
    dpKdT(23) = 0.006666667d0
    dpKdT(24) = -0.008333333d0
    dpKdT(25) = -0.055833333d0
    dpKdT(26) = -0.031932573d0 ! water

    end subroutine read_all

    end module read_parameters
