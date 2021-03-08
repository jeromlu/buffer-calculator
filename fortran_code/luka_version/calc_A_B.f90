    module calc_A_B
    ! this module calculates the values of A and B buffers, as well as the mixtures

    ! use
    use read_parameters
    use calc_grams

    ! parameters
    implicit none

    ! subroutines
    contains

    subroutine calc_buffs
    ! this gets the data and sends it to the workers!

    ! parameters
    implicit none
    double precision :: molarities(11), pH ,grams(20,2), conduc(2), elec_error(2), ion_str(2), step_size, perc_B
    double precision :: mol_mix(11), Na_mix, Cl_mix, cond_mix, error_mix, Na_pres(2), Cl_pres(2)
    double precision :: temp_min, temp_step, temp_mid, temp_test, visc_temp
    double precision, allocatable :: ion_mix(:), pH_mix(:)
    integer :: phos_type, loop_calc,loop_grams, n_steps, loop_mix, loop_temp, loop_diff

    ! calculations
    ! opening interesting files
    open(unit=25, file='out_A&B_sngl.dat', action='write', status='replace') ! A and B related values!

    ! getting A and B gram values!
    do loop_calc = 1, 2

        ! setting input the values
        molarities(:) = input_molarities(:,loop_calc) ! setting the proper molarities
        pH = input_pH(loop_calc) ! getting the proper pH
        phos_type = input_phos_type(loop_calc) ! setting the hosphate type

        ! Calling the routine.
        call calc_all_grams(molarities, pH, phos_type) ! getting the pH values!

        ! Setting output values!
        grams(:,loop_calc) = out_grams(:)* input_V(loop_calc)
        ion_str(loop_calc) = ion_str_real 
        conduc(loop_calc) = conductivity
        elec_error(loop_calc) = electro_error**2d0
        Na_pres(loop_calc) = Na_needed
        Cl_pres(loop_calc) = Cl_needed
    end do

    ! Write to file 'out_A&B_sngl.dat'
    do loop_grams = 1,20
        write(25,5) grams (loop_grams,:) ! Writing gram values.
    end do
    write(25,5) ion_str
    write(25,5) conduc
    write(25,5) elec_error
    
    ! getting everything in between!
    perc_B = 0d0
    n_steps = 30! number of steps to try! 
    step_size = 100d0 / dble(n_steps-1)
    perc_B = perc_B
    allocate (pH_mix(n_steps), ion_mix(n_steps))
    pH = input_pH(1)
    
    open(unit=26, file='out_mix.dat', action='write', status='replace') ! mix related values!
    
    do loop_mix = 0,n_steps-1
        perc_B = step_size * dble (loop_mix) ! getting the percent of B
        ! getting molarities
        mol_mix(:) = (1d0 - perc_B/100d0) * input_molarities(:,1) + (perc_B/100d0) * input_molarities(:,2)
        ! and Na and Cl values
        Na_mix = (1d0 - perc_B/100d0) * Na_pres(1) + (perc_B/100d0) * Na_pres(2)
        Cl_mix = (1d0 - perc_B/100d0) * Cl_pres(1) + (perc_B/100d0) * Cl_pres(2)
        call calc_all_pH (mol_mix, pH, Na_mix, Cl_mix)
        pH_mix(loop_mix+1) = pH_out
        ion_mix(loop_mix+1) = ion_str_real 
        cond_mix = conductivity
        error_mix = electro_error**2d0
        ! writing results        
        write(26,5) perc_B, ion_mix(loop_mix+1), pH_mix(loop_mix+1), cond_mix, error_mix
        pH = pH_mix(loop_mix+1)
    end do
    
    ! getting temperature dependency
    n_steps = 10 ! number of steps to try
    ! temperature values
    Temp_min = 298d0 - 10d0 ! in K  
    Temp_step = 2d0 ! in K
    Temp_mid = 298d0 ! in K
    ! opening files 
    open(unit=55, file='out_temp.dat', action='write', status = 'replace') ! temperature values! 
    
    
    do loop_temp = 0,n_steps
        Temp_test = temp_min + temp_step*(loop_temp) ! setting temperatures for new calculations
        ! updating the pKa values
        pKa_acet = 4.76d0 + dpKdT(1)*(-Temp_mid+Temp_test)

        pKa_phos(1) = 2.15d0+ dpKdT(2)*(-Temp_mid+Temp_test)
        pKa_phos(2) = 7.21d0+ dpKdT(3)*(-Temp_mid+Temp_test)
        pKa_phos(3) = 12.33d0+ dpKdT(4)*(-Temp_mid+Temp_test)

        pKa_mala(1) = 3.4d0+ dpKdT(5)*(-Temp_mid+Temp_test)
        pKa_mala(2) = 5.13d0+ dpKdT(6)*(-Temp_mid+Temp_test)

        pKa_succ(1) = 4.21d0+ dpKdT(7)*(-Temp_mid+Temp_test)
        pKa_succ(2) = 5.64d0+ dpKdT(8)*(-Temp_mid+Temp_test)

        pKa_tris = 8.07d0+ dpKdT(9)*(-Temp_mid+Temp_test)

        pKa_citr(1) = 3.13d0+ dpKdT(10)*(-Temp_mid+Temp_test)
        pKa_citr(2) = 4.76d0+ dpKdT(11)*(-Temp_mid+Temp_test)
        pKa_citr(3) = 6.4d0+ dpKdT(12)*(-Temp_mid+Temp_test)

        pKa_argi(1) = 2.09d0+ dpKdT(13)*(-Temp_mid+Temp_test)
        pKa_argi(2) = 9d0+ dpKdT(14)*(-Temp_mid+Temp_test)
        pKa_argi(3) = 12.1d0+ dpKdT(15)*(-Temp_mid+Temp_test)

        pKa_hist(1) = 1.7d0+ dpKdT(16)*(-Temp_mid+Temp_test)
        pKa_hist(2) = 9.09d0+ dpKdT(17)*(-Temp_mid+Temp_test)
        pKa_hist(3) = 6.04d0+ dpKdT(18)*(-Temp_mid+Temp_test)
    
        pKa_cyst(1) = 1.91d0+ dpKdT(19)*(-Temp_mid+Temp_test)
        pKa_cyst(2) = 10.28d0+ dpKdT(20)*(-Temp_mid+Temp_test)
        pKa_cyst(3) = 8.14d0+ dpKdT(21)*(-Temp_mid+Temp_test)

        pKa_EDTA(1) = 2d0+ dpKdT(22)*(-Temp_mid+Temp_test)
        pKa_EDTA(2) = 2.7d0+ dpKdT(23)*(-Temp_mid+Temp_test)
        pKa_EDTA(3) = 6.2d0+ dpKdT(24)*(-Temp_mid+Temp_test)
        pKa_EDTA(4) = 10.31d0+ dpKdT(25)*(-Temp_mid+Temp_test)
        
        pKw = 14d0+ dpKdT(26)*(-Temp_mid+Temp_test)
        
        ! updating the debye huckel paramter A values
        A_DH = A_DH_relationship(1)*Temp_test + A_DH_relationship(2)
        
        ! updating the diffusivities
        ! viscosity
        visc_temp = visc(1) * Temp_test + visc(2)
        do loop_diff = 1,23
            Diff(loop_diff) = kB6pi * Temp_test/ion_rad(loop_diff)/visc_temp
        end do
       ! write(*,*) temp_test
       ! write(*,*) Diff
        !pause
        ! calling the pH solver! 
        call calc_all_pH (input_molarities(:,1), input_pH(1), Na_pres(1), Cl_pres(1))
        
        ! writting results
        write(55,5) Temp_test, ion_str_real, pH_out, conductivity, electro_error**2d0
        
        
    end do
    

    close (25)
    close (26)
    close (55)
5   format(15(f30.5))

    end subroutine calc_buffs
    end module calc_A_B



