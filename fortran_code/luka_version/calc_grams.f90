    module calc_grams
    ! given a pH value and molarities, will calculate the amounts of the different components needed!

    ! use stuff
    use read_parameters

    ! parameters
    implicit none
    double precision :: ion_str_error, Na_needed, Cl_needed, electro_error, ion_str_avail, ion_str_real
    double precision :: AcO, H2PO4, HPO4, PO4, succH, succ, tris, H2cit, Hcit, cit, malaH, mala, arg1, arg2, arg3
    double precision :: his1, his2, his3, cys1, cys2, cys3, H3EDTA, H2EDTA, HEDTA, EDTA, H, OH
    double precision :: out_grams(20), conductivity
    double precision :: in_molarities(11) , in_pH, pH_out
    integer :: in_phos_type, loop_in
    integer :: ion_break , break

    ! subroutines
    contains

    subroutine calc_all_pH (molarities, pH, Na, Cl)
    ! gets the pH values of the mixtures between A and B, as well as conductivities, ionic strengths, ...

    ! parameters
    implicit none
    double precision :: molarities(11), pH,  Na, Cl, ion_guess

    ! calculations
    in_molarities (:) = molarities(:)

    ! getting the pH value!
    pH_out = get_pH ( pH, Na, Cl)

    ! getting the ionic strength!
    ion_guess = 0.5d0 * (Na+Cl)
    ion_str_real = get_ion_str_real (ion_guess, Na, Cl)

    ! getting the conductivity
    call calc_conductivity (Na, Cl)

    end subroutine calc_all_pH

    double precision function get_pH ( pH_guess, Na_val, Cl_val)
    ! calculates the pH!

    ! parameters
    implicit none
    double precision pH_guess, Na_val, Cl_val, error_pH, error_pH2, pH_guess2, ion_guess, slope
    integer :: counter, break_pH

    ! calculations
    ! initializing
    counter = 1
    break_pH = 0

    ! we first calculate the ionic strength
    ion_guess = 0.5d0 * (Na_val + Cl_val)

    ! getting the initial error
    in_pH = pH_guess
    call calc_electroneutrality_error(ion_guess, Na_val, Cl_val)
    error_pH = electro_error

    if (abs(error_pH) .lt. 1d-15) then
        break_pH = 1 ! we have convergence!
        get_pH = pH_guess
    else
        ! getting the other initial error
        if (error_pH .gt. 0d0) then ! in this case we have too much H+
            pH_guess2 = pH_guess - (0.1d0)
        else ! in this case we have not enough H+
            pH_guess2 = pH_guess + (0.1d0)
        end if

        in_pH = pH_guess2
        call calc_electroneutrality_error(ion_guess, Na_val, Cl_val)
        error_pH2 = electro_error

        ! we use a derivative method here
        do while( break_pH .eq. 0)
            slope = (error_pH2 - error_pH) / (pH_guess2 - pH_guess)
            pH_guess = pH_guess2
            error_pH = error_pH2
            pH_guess2 = -error_pH2/slope + pH_guess
            in_pH = pH_guess2
            call calc_electroneutrality_error(ion_guess, Na_val, Cl_val)
            error_pH2 = electro_error

            ! checking for convergence
            if (abs(error_pH2) .lt. 1d-15) then
                break_pH = 1 ! we have convergence!
            end if
            if (counter .gt. 1000) then
                break_pH = -1 ! no convergence
            end if
            counter = counter +1
        end do
        get_pH = pH_guess2
    end if
    end function get_pH

    subroutine calc_all_grams(molarities, pH, phos_type )
    ! calculates all in this sheet

    ! parameters
    implicit none
    double precision :: molarities(11), pH
    integer :: phos_type

    ! setting the values
    in_molarities (:) = molarities(:)
    in_pH = pH
    in_phos_type = phos_type

    ! now we need to assign the different amounts of Na and Cl to the different species!
    call calc_Na_Cl_needed

    ! next step is to assign the required valueS!
    call assign_grams

    ! getting the real ionic strength!
    ion_str_real = get_ion_str_real (0.5d0 * Na_needed + 0.5d0 * Cl_needed, Na_needed, Cl_needed)

    ! then we get the conductivities
    call calc_conductivity (Na_needed, Cl_needed)

    end subroutine calc_all_grams

    subroutine calc_conductivity( Na_in, Cl_in)
    ! this subroutine calculates the conductivity!

    ! parameters
    implicit none
    double precision :: cond(23), gamma , alpha,z, Na_in, Cl_in
    integer :: counter

    ! first we get
    counter = 1
    ! for AcO-
    z=-1d0
    gamma = 10d0**(-(z)**2d0*(A_DH*sqrt(ion_str_real)/(1+sqrt(ion_str_real))-B_DH(1)*ion_str_real))
    if (ion_str_real .lt. 0.36d0*abs(z)) then
        alpha = 0.6d0/abs(z)
    else
        alpha = sqrt( ion_str_real/abs(z))
    end if
    cond(counter) = Diff(counter)*(z)**2d0*gamma**alpha*AcO

    ! for H2PO4-
    counter = counter + 1
    z=-1d0
    gamma = 10d0**(-(z)**2d0*(A_DH*sqrt(ion_str_real)/(1+sqrt(ion_str_real))-B_DH(2)*ion_str_real))
    if (ion_str_real .lt. 0.36d0*abs(z)) then
        alpha = 0.6d0/abs(z)
    else
        alpha = sqrt( ion_str_real/abs(z))
    end if
    cond(counter) = Diff(counter)*(z)**2d0*gamma**alpha*H2PO4

    ! for HPO4 2-
    counter = counter + 1
    z=-2d0
    gamma = 10d0**(-(z)**2d0*(A_DH*sqrt(ion_str_real)/(1+sqrt(ion_str_real))-B_DH(2)*ion_str_real))
    if (ion_str_real .lt. 0.36d0*abs(z)) then
        alpha = 0.6d0/abs(z)
    else
        alpha = sqrt( ion_str_real/abs(z))
    end if
    cond(counter) = Diff(counter)*(z)**2d0*gamma**alpha*HPO4

    ! for PO4 3-
    counter = counter + 1
    z=-3d0
    gamma = 10d0**(-(z)**2d0*(A_DH*sqrt(ion_str_real)/(1+sqrt(ion_str_real))-B_DH(2)*ion_str_real))
    if (ion_str_real .lt. 0.36d0*abs(z)) then
        alpha = 0.6d0/abs(z)
    else
        alpha = sqrt( ion_str_real/abs(z))
    end if
    cond(counter) = Diff(counter)*(z)**2d0*gamma**alpha*PO4

    ! for succH -
    counter = counter + 1
    z=-1d0
    gamma = 10d0**(-(z)**2d0*(A_DH*sqrt(ion_str_real)/(1+sqrt(ion_str_real))-B_DH(3)*ion_str_real))
    if (ion_str_real .lt. 0.36d0*abs(z)) then
        alpha = 0.6d0/abs(z)
    else
        alpha = sqrt( ion_str_real/abs(z))
    end if
    cond(counter) = Diff(counter)*(z)**2d0*gamma**alpha*succH

    ! for succ 2 -
    counter = counter + 1
    z=-2d0
    gamma = 10d0**(-(z)**2d0*(A_DH*sqrt(ion_str_real)/(1+sqrt(ion_str_real))-B_DH(3)*ion_str_real))
    if (ion_str_real .lt. 0.36d0*abs(z)) then
        alpha = 0.6d0/abs(z)
    else
        alpha = sqrt( ion_str_real/abs(z))
    end if
    cond(counter) = Diff(counter)*(z)**2d0*gamma**alpha*succ

    ! for tris -
    counter = counter + 1
    z=1d0
    gamma = 10d0**(-(z)**2d0*(A_DH*sqrt(ion_str_real)/(1+sqrt(ion_str_real))-B_DH(4)*ion_str_real))
    if (ion_str_real .lt. 0.36d0*abs(z)) then
        alpha = 0.6d0/abs(z)
    else
        alpha = sqrt( ion_str_real/abs(z))
    end if
    cond(counter) = Diff(counter)*(z)**2d0*gamma**alpha*tris

    ! for H2cit -
    counter = counter + 1
    z=-1d0
    gamma = 10d0**(-(z)**2d0*(A_DH*sqrt(ion_str_real)/(1+sqrt(ion_str_real))-B_DH(5)*ion_str_real))
    if (ion_str_real .lt. 0.36d0*abs(z)) then
        alpha = 0.6d0/abs(z)
    else
        alpha = sqrt( ion_str_real/abs(z))
    end if
    cond(counter) = Diff(counter)*(z)**2d0*gamma**alpha*H2cit

    ! for Hcit2 -
    counter = counter + 1
    z=-2d0
    gamma = 10d0**(-(z)**2d0*(A_DH*sqrt(ion_str_real)/(1+sqrt(ion_str_real))-B_DH(5)*ion_str_real))
    if (ion_str_real .lt. 0.36d0*abs(z)) then
        alpha = 0.6d0/abs(z)
    else
        alpha = sqrt( ion_str_real/abs(z))
    end if
    cond(counter) = Diff(counter)*(z)**2d0*gamma**alpha*Hcit

    ! for cit 3-
    counter = counter + 1
    z=-3d0
    gamma = 10d0**(-(z)**2d0*(A_DH*sqrt(ion_str_real)/(1+sqrt(ion_str_real))-B_DH(5)*ion_str_real))
    if (ion_str_real .lt. 0.36d0*abs(z)) then
        alpha = 0.6d0/abs(z)
    else
        alpha = sqrt( ion_str_real/abs(z))
    end if
    cond(counter) = Diff(counter)*(z)**2d0*gamma**alpha*cit

    ! for malaH
    counter = counter + 1
    z=-1d0
    gamma = 10d0**(-(z)**2d0*(A_DH*sqrt(ion_str_real)/(1+sqrt(ion_str_real))-B_DH(6)*ion_str_real))
    if (ion_str_real .lt. 0.36d0*abs(z)) then
        alpha = 0.6d0/abs(z)
    else
        alpha = sqrt( ion_str_real/abs(z))
    end if
    cond(counter) = Diff(counter)*(z)**2d0*gamma**alpha*malah

    ! for mala
    counter = counter + 1
    z=-2d0
    gamma = 10d0**(-(z)**2d0*(A_DH*sqrt(ion_str_real)/(1+sqrt(ion_str_real))-B_DH(6)*ion_str_real))
    if (ion_str_real .lt. 0.36d0*abs(z)) then
        alpha = 0.6d0/abs(z)
    else
        alpha = sqrt( ion_str_real/abs(z))
    end if
    cond(counter) = Diff(counter)*(z)**2d0*gamma**alpha*mala

    ! for arg
    counter = counter + 1
    z=+1d0
    gamma = 10d0**(-(z)**2d0*(A_DH*sqrt(ion_str_real)/(1+sqrt(ion_str_real))-B_DH(7)*ion_str_real))
    if (ion_str_real .lt. 0.36d0*abs(z)) then
        alpha = 0.6d0/abs(z)
    else
        alpha = sqrt( ion_str_real/abs(z))
    end if
    cond(counter) = Diff(counter)*(z)**2d0*gamma**alpha*arg3

    ! for his
    counter = counter + 1
    z=+1d0
    gamma = 10d0**(-(z)**2d0*(A_DH*sqrt(ion_str_real)/(1+sqrt(ion_str_real))-B_DH(8)*ion_str_real))
    if (ion_str_real .lt. 0.36d0*abs(z)) then
        alpha = 0.6d0/abs(z)
    else
        alpha = sqrt( ion_str_real/abs(z))
    end if
    cond(counter) = Diff(counter)*(z)**2d0*gamma**alpha*his3

    ! for cys
    counter = counter + 1
    z=-1d0
    gamma = 10d0**(-(z)**2d0*(A_DH*sqrt(ion_str_real)/(1+sqrt(ion_str_real))-B_DH(9)*ion_str_real))
    if (ion_str_real .lt. 0.36d0*abs(z)) then
        alpha = 0.6d0/abs(z)
    else
        alpha = sqrt( ion_str_real/abs(z))
    end if
    cond(counter) = Diff(counter)*(z)**2d0*gamma**alpha*cys3

    ! for H3EDTA
    counter = counter + 1
    z=-1d0
    gamma = 10d0**(-(z)**2d0*(A_DH*sqrt(ion_str_real)/(1+sqrt(ion_str_real))-B_DH(10)*ion_str_real))
    if (ion_str_real .lt. 0.36d0*abs(z)) then
        alpha = 0.6d0/abs(z)
    else
        alpha = sqrt( ion_str_real/abs(z))
    end if
    cond(counter) = Diff(counter)*(z)**2d0*gamma**alpha*H3EDTA

    ! for H2 EDTA 2-
    counter = counter + 1
    z=-2d0
    gamma = 10d0**(-(z)**2d0*(A_DH*sqrt(ion_str_real)/(1+sqrt(ion_str_real))-B_DH(10)*ion_str_real))
    if (ion_str_real .lt. 0.36d0*abs(z)) then
        alpha = 0.6d0/abs(z)
    else
        alpha = sqrt( ion_str_real/abs(z))
    end if
    cond(counter) = Diff(counter)*(z)**2d0*gamma**alpha*H2EDTA

    ! for H EDTA 3-
    counter = counter + 1
    z=-3d0
    gamma = 10d0**(-(z)**2d0*(A_DH*sqrt(ion_str_real)/(1+sqrt(ion_str_real))-B_DH(10)*ion_str_real))
    if (ion_str_real .lt. 0.36d0*abs(z)) then
        alpha = 0.6d0/abs(z)
    else
        alpha = sqrt( ion_str_real/abs(z))
    end if
    cond(counter) = Diff(counter)*(z)**2d0*gamma**alpha*HEDTA

    ! for  EDTA 4-
    counter = counter + 1
    z=-4d0
    gamma = 10d0**(-(z)**2d0*(A_DH*sqrt(ion_str_real)/(1+sqrt(ion_str_real))-B_DH(10)*ion_str_real))
    if (ion_str_real .lt. 0.36d0*abs(z)) then
        alpha = 0.6d0/abs(z)
    else
        alpha = sqrt( ion_str_real/abs(z))
    end if
    cond(counter) = Diff(counter)*(z)**2d0*gamma**alpha*EDTA

    ! for  Na+
    counter = counter + 1
    z=+1d0
    gamma = 10d0**(-(z)**2d0*(A_DH*sqrt(ion_str_real)/(1+sqrt(ion_str_real))-B_DH(10)*ion_str_real))
    if (ion_str_real .lt. 0.36d0*abs(z)) then
        alpha = 0.6d0/abs(z)
    else
        alpha = sqrt( ion_str_real/abs(z))
    end if
    cond(counter) = Diff(counter)*(z)**2d0*gamma**alpha*Na_in

    ! for  CL-
    counter = counter + 1
    z=-1d0
    gamma = 10d0**(-(z)**2d0*(A_DH*sqrt(ion_str_real)/(1+sqrt(ion_str_real))-B_DH(10)*ion_str_real))
    if (ion_str_real .lt. 0.36d0*abs(z)) then
        alpha = 0.6d0/abs(z)
    else
        alpha = sqrt( ion_str_real/abs(z))
    end if
    cond(counter) = Diff(counter)*(z)**2d0*gamma**alpha*Cl_in

    ! forH+
    counter = counter + 1
    z=1d0
    gamma = 10d0**(-(z)**2d0*(A_DH*sqrt(ion_str_real)/(1+sqrt(ion_str_real))-B_DH(10)*ion_str_real))
    if (ion_str_real .lt. 0.36d0*abs(z)) then
        alpha = 0.6d0/abs(z)
    else
        alpha = sqrt( ion_str_real/abs(z))
    end if
    cond(counter) = Diff(counter)*(z)**2d0*gamma**alpha*10**-in_pH

    ! for OH-
    counter = counter + 1
    z=1d0
    gamma = 10d0**(-(z)**2d0*(A_DH*sqrt(ion_str_real)/(1+sqrt(ion_str_real))-B_DH(10)*ion_str_real))
    if (ion_str_real .lt. 0.36d0*abs(z)) then
        alpha = 0.6d0/abs(z)
    else
        alpha = sqrt( ion_str_real/abs(z))
    end if
    cond(counter) = Diff(counter)*(z)**2d0*gamma**alpha*((10**-pKw)/10**-in_pH)

    ! and finally the full conductivity!
    conductivity = 3755400 * sum(cond)*10d0*1000d0

    end subroutine calc_conductivity

    subroutine assign_grams
    ! this subroutine takes the value of Na_needed and Cl_needed and assigns amounts of grams og the buffering components, as well as the amounts of Naoh and Hcl !

    ! paramters
    implicit none
    double precision :: extra_Na, Na_counter, Na_ac_max, Na_phos_max, Na_ac_min, Na_phos_min, Na_ac_alloc, Na_phos_alloc
    integer :: Na_ac_allowed, Na_phos_allowed, cl_tris_allowed, cl_hist_allowed
    double precision :: extra_cl, Cl_counter, cl_tris_max,cl_hist_max, tris_Hcl, tris_base, hist_hcl, hist_base
    double precision :: NaAcO, AcOH, H3PO4, NaH2PO4, Na2HPO4, Na3PO4

    ! calculating!
    ! known species (these are the species that are inputed without any Na or Cl(i.e. succinate, citric acid....)
    ! checking for availability of these species
    out_grams(8) = in_molarities(3) * MW(7) ! succinate
    out_grams(11) = in_molarities(5) * MW(17)  ! citric acid
    out_grams(12) = in_molarities(6) * MW(18)! malicacid
    out_grams(13) = in_molarities(7) * MW(10)! arginine HcL
    out_grams(16) = in_molarities(9) * MW(19)! cysteine
    out_grams(17) = in_molarities(10) * MW(20)! EDTA
    out_grams(18) = in_molarities(11) * MW(14)! NaCl

    ! Na +
    extra_Na = Na_needed - in_molarities (11) - 2d0*in_molarities(10) ! this is the amount of Na+ not provided for by NaCl or EDTA

    ! two species contribute to the  Na+ balance : Acetate (NaAcO-) and phsophate (depends on which species it is, i.e. NaH2PO4)
    ! chekcing for the availability of these
    Na_counter = 0d0
    if (in_molarities(1) .gt. 0d0) then
        Na_counter = Na_counter + 1d0
        Na_ac_allowed = 1
    else
        Na_ac_allowed = 0
    end if

    if (in_molarities(2) .gt. 0d0) then
        Na_counter = Na_counter + 1d0
        Na_phos_allowed = 1
    else
        Na_phos_allowed = 0
    end if

    if (Na_counter .eq. 0d0) then
        NaAcO = 0d0
        AcOH = 0d0
        NaH2PO4 = 0d0
        H3PO4  = 0d0
        Na2HPO4 = 0d0
        Na3PO4 = 0d0

    else
        ! max Na provided by acetate & phos
        Na_ac_max = in_molarities(1) ! if acetate is allowed, then max acetate is full NaAcO
        Na_phos_max = in_molarities(2) * dble(in_phos_type )! if phosphate is allowed, then we full full phos_type (NaH2, Na2H or Na3) PO4
        if (Na_counter .eq. 1d0) then ! only one species allowed
            if (Na_ac_allowed .eq. 1d0) then ! only acetate allowed

                if (extra_Na .gt. Na_ac_max) then ! more Na needed that NaAc is available
                    NaAcO = Na_ac_max
                    AcOH = 0d0 ! we don't have this since we need all the actate Na we can get!
                else
                    NaAcO = extra_Na ! we only need the described Na value
                    AcOH = in_molarities(1) - NaAcO
                end if
            end if

            if (Na_phos_allowed .eq. 1d0) then ! only phosphjate is allower
                if (extra_Na .gt. Na_phos_max) then ! more Na needed that (Na)n-phos is available
                    if (in_phos_type  .eq. 1) then
                        NaH2PO4 = Na_phos_max / dble(in_phos_type)
                        H3PO4  = 0d0
                        Na2HPO4 = 0d0
                        Na3PO4 =0d0
                    else if (in_phos_type .eq. 2) then
                        Na2HPO4 = Na_phos_max / dble(in_phos_type)
                        H3PO4  = 0d0
                        NaH2PO4 = 0d0
                        Na3PO4 =0d0
                    else if (in_phos_type .eq. 3) then
                        Na3PO4 = Na_phos_max / dble(in_phos_type)
                        H3PO4  = 0d0
                        NaH2PO4 = 0d0
                        Na2HPO4 =0d0
                    end if
                else ! we dont need all the phos to be Na max type
                    if (in_phos_type  .eq. 1) then
                        NaH2PO4 = extra_Na ! we only need the described Na value
                        H3PO4  = in_molarities(2) - NaH2PO4
                        Na2HPO4 = 0d0
                        Na3PO4 =0d0
                    else if (in_phos_type .eq. 2) then
                        Na2HPO4 = extra_Na - in_molarities(2)
                        H3PO4  = 0d0
                        NaH2PO4 = in_molarities(2) - Na2HPO4
                        Na3PO4 =0d0
                    else if (in_phos_type .eq. 3) then
                        Na3PO4 = extra_Na - 2d0 * in_molarities(2)
                        H3PO4  = 0d0
                        NaH2PO4 = 0d0
                        Na2HPO4 = in_molarities(2) - Na3PO4
                    end if
                end if
            end if

        else ! (Na_counter must be 2d0! and both acetate and phosphate are allowed
            if (extra_Na .gt. Na_ac_max + Na_phos_max) then ! we need all the Na we can get!
                if (in_phos_type  .eq. 1) then
                    NaH2PO4 = Na_phos_max / dble(in_phos_type)
                    H3PO4  = 0d0
                    Na2HPO4 = 0d0
                    Na3PO4 =0d0
                else if (in_phos_type .eq. 2) then
                    Na2HPO4 = Na_phos_max / dble(in_phos_type)
                    H3PO4  = 0d0
                    NaH2PO4 = 0d0
                    Na3PO4 =0d0
                else if (in_phos_type .eq. 3) then
                    Na3PO4 = Na_phos_max / dble(in_phos_type)
                    H3PO4  = 0d0
                    NaH2PO4 = 0d0
                    Na2HPO4 =0d0
                end if
                NaAcO = Na_ac_max
                AcOH  = 0d0
            end if

            if(extra_Na .le. Na_ac_max + Na_phos_max) then ! we have enough Na with phos and ac alone!
                if (in_phos_type  .eq. 1) then
                    NaH2PO4 = Na_phos_max / dble(in_phos_type)
                    H3PO4  = 0d0
                    Na2HPO4 = 0d0
                    Na3PO4 =0d0
                else if (in_phos_type .eq. 2) then
                    Na2HPO4 = Na_phos_max / dble(in_phos_type)
                    H3PO4  = 0d0
                    NaH2PO4 = 0d0
                    Na3PO4 =0d0
                else if (in_phos_type .eq. 3) then
                    Na3PO4 = Na_phos_max / dble(in_phos_type)
                    H3PO4  = 0d0
                    NaH2PO4 = 0d0
                    Na2HPO4 =0d0
                end if
                NaAcO = extra_Na - NaH2PO4 - 2d0 * Na2HPO4 - 3d0 * Na3PO4
                AcOH  = in_molarities(1) - NaAcO
            end if

            if (extra_Na .le. Na_phos_max) then ! we have enough Na in phosphate!
                if (in_phos_type  .eq. 1) then
                    NaH2PO4 = extra_Na ! we only need the described Na value
                    H3PO4  = in_molarities(2) - NaH2PO4
                    Na2HPO4 = 0d0
                    Na3PO4 =0d0
                else if (in_phos_type .eq. 2) then
                    Na2HPO4 = extra_Na - in_molarities(2)
                    H3PO4  = 0d0
                    NaH2PO4 = in_molarities(2) - Na2HPO4
                    Na3PO4 =0d0
                else if (in_phos_type .eq. 3) then
                    Na3PO4 = extra_Na - 2d0 * in_molarities(2)
                    H3PO4  = 0d0
                    NaH2PO4 = 0d0
                    Na2HPO4 = in_molarities(2) - Na3PO4
                end if
                NaAcO = 0d0
                AcOH  = in_molarities(1)
            end if

        end if
    end if
















    ! finally we can get NaOH is needed
    out_grams(19) = (extra_Na - NaAcO - NaH2PO4 - 2*Na2HPO4 - 3* Na3PO4) /10d0 * dens(3) ! mol/L needed /10 mol/L *g/L * 1L

    ! and we can also get the amounts of the acetate and phosphate needed
    out_grams(3) = NaAcO * MW(2) ! NaAcO
    out_grams(2) = AcOH * MW(1) ! AcOH
    out_grams(4) = H3PO4 * MW(3)/0.85d0 ! H3PO4
    out_grams(5) = NaH2PO4 * MW(4) ! NaH2PO4
    out_grams(6) = NA2HPO4 * MW(5) ! Na2HPO4
    out_grams(7) = NA3PO4 * MW(6) ! NA3PO4

    ! for Cl need also to consider cystein HcL monohydrate
    extra_Cl = Cl_needed - (in_molarities(11) + in_molarities(7) + in_molarities(9)) ! removing NaCl, arginine, and cysteine from required Cl amounts

    ! two species contribute to the  Cl+ balance : tris and histidine
    ! chekcing for the availability of these
    Cl_counter = 0d0
    if (in_molarities(4) .gt. 0d0) then
        Cl_counter = Cl_counter + 1d0
        Cl_tris_allowed = 1
    else
        Cl_tris_allowed = 0
    end if
    if (in_molarities(8) .gt. 0d0) then
        Cl_counter = Cl_counter + 1d0
        Cl_hist_allowed = 1
    else
        Cl_hist_allowed = 0
    end if

    ! getting minimum and maximum tris and istidine amounts!
    Cl_tris_max = in_molarities(4) ! if all is tris HCl
    Cl_hist_max = in_molarities(8) ! if all his HCl

    ! assigning the amounts
    if (Cl_counter .eq. 0d0) then
        tris_HCl = 0d0
        tris_base = 0d0
        hist_Hcl = 0d0
        hist_base = 0d0
    else
        if (extra_Cl / Cl_counter .gt. Cl_tris_max) then
            tris_HCl = Cl_tris_max
        else
            tris_HCl = extra_Cl / Cl_counter
        end if
        tris_base = in_molarities(4) - tris_HCl
        if (extra_Cl / Cl_counter .gt. Cl_hist_max) then
            hist_HCl = Cl_hist_max
        else
            Hist_HCl = extra_Cl / Cl_counter
        end if
        hist_base = in_molarities(8) - hist_HCl
    end if

    ! and getting the amount of HCL, if needed
    out_grams(20) = (extra_Cl - tris_HCl - hist_HCL) * MW(16)/.25d0

    ! getting the grams...
    out_grams(9) = tris_base *MW(8) ! tris base
    out_grams(10) = tris_HCl * MW (9) ! tris HCL
    out_grams(14) = hist_base * MW(11) ! histidine base
    out_grams(15) = hist_HCl * MW(12) ! histidine HCL

    ! and finally for water!
    out_grams(1) = (1-(out_grams(2)/dens(2) + out_grams(4)/dens(5) + out_grams(5)*MW(15)/MW(4)/dens(1) &
        + out_grams(15)*MW(15)/MW(12)/dens(1) + out_grams(16)*MW(15)/MW(19)/dens(1)+ out_grams(19)/dens(3) &
        + out_grams(20)/dens(4)+2d0*out_grams(17)*MW(15)/MW(20)/dens(1)+out_grams(11)*MW(15)/MW(17)/dens(1))) *dens(1)

    end subroutine assign_grams

    subroutine calc_Na_Cl_needed
    ! gets the required amounts of Na and Cl

    ! parmateters
    implicit none
    double precision :: Na_min , Cl_min
    double precision :: ion_str

    ! getting the required Na and Cl values!
    ! the minimum value of Na is
    Na_min = in_molarities(11) + (dble(in_phos_type)-1d0)*in_molarities(2) + 2d0*in_molarities(10) ! NaCl + if we have phophate buffer we need to consider we might have some of it in there and also EDTA!!
    Cl_min = in_molarities(11) + in_molarities(7) + in_molarities(9) ! NaCl and adding any of arginine or cystein will automatically bring in some extra Cl!

    ! the minimum ionic strength is then
    ion_str = 0.5*(Na_min + Cl_min)

    ! getting the electroneutrality at the initial values
    call calc_electroneutrality_error(ion_str, Na_min, Cl_min)
    
    if (electro_error .lt. 0d0) then ! we need Na+ in this case!
        Cl_needed = Cl_min
        Na_needed = get_Na_conc(Na_min)
    else if (electro_error .eq. 0d0) then ! we dont need Na or Cl!
        Na_needed = Na_min
        Cl_needed = Cl_min
    else if (electro_error .gt. 0d0 ) then ! we need Cl- in this case!
        Na_needed = Na_min
        Cl_needed = get_Cl_conc(Cl_min)
    end if

    ! we then want to get the final values of ionic strength and of the individual valueS!
    ion_str_avail = 0.5d0 * (Na_Needed + Cl_needed)
    call calc_species_conc(ion_str_avail)

    end subroutine calc_Na_Cl_needed

    double precision function get_Na_conc(Na_guess)
    ! gets the value of Na that is required!

    !parameters
    implicit none
    double precision :: Na_guess, ion_guess
    double precision :: Na_guess2, error_Na, error_Na2, slope
    integer :: counter

    ! initialising
    counter = 1
    break = 0
    error_Na = electro_error
    ! setting new Na_guess value
    Na_guess2 = 2d0*Na_guess + 30d-3
    ! we can then get the ionic strength for these values of Na and Cl
    ion_guess = 0.5d0 * Na_guess2 + 0.5d0 * Cl_needed
    ! and call the electroneutrality equation again
    call calc_electroneutrality_error(ion_guess, Na_guess2, Cl_needed)
    error_Na2 = electro_error
    ! we then use a derivative method
    do while (break .eq. 0)
        slope = (error_Na2-error_Na) / (Na_guess2 - Na_guess)
        Na_guess = Na_guess2
        error_Na = error_Na2
        Na_guess2 = -error_Na2/slope + Na_guess
        ! checking for neg values
        if (Na_guess2 .lt. 0d0) then
            Na_guess2 = 0d0     
        end if        
        ! we can then get the ionic strength for these values of Na and Cl
        ion_guess = 0.5d0 * Na_guess2 + 0.5d0 * Cl_needed
        ! and call the elctroneutrality again
        call calc_electroneutrality_error(ion_guess, Na_guess2, Cl_needed)
        error_Na2 = electro_error

        ! checking for convergence
        if (error_Na2 .eq. 1d-25) then
            break = 1 ! we have convergence!
        end if
        if (abs(Na_guess/Na_guess2 -1d0) .lt. 1d-25) then
            break = 2 ! no more change in the value of x!
        end if
        if (abs(error_Na/error_Na2 - 1d0) .lt. 1d-25) then
            break = 3 ! no more change in the value of f(x)
        end if
        if (counter .gt. 1000) then
            break = -1 ! no convergence
        end if
        counter = counter +1
    end do
    get_Na_conc = Na_guess2

    end function get_Na_conc

    double precision function get_Cl_conc(Cl_guess)
    ! gets the value of Na that is required!

    !parameters
    implicit none
    double precision :: Cl_guess, ion_guess
    double precision :: Cl_guess2, error_Cl, error_Cl2, slope
    integer :: counter

    ! initialising
    counter = 1
    break = 0
    error_Cl = electro_error
    ! setting new Cl_guess value
    Cl_guess2 = 2d0*Cl_guess + 30d-3
    ! we can then get the ionic strength for these values of Na and Cl
    ion_guess = 0.5d0 * Na_needed + 0.5d0 * Cl_guess2
    ! and call the electroneutrality equation again
    call calc_electroneutrality_error(ion_guess, Na_needed, Cl_guess2)
    error_Cl2 = electro_error
    ! we then use a derivative method
    do while (break .eq. 0)
        slope = (error_Cl2-error_Cl) / (Cl_guess2 - Cl_guess)
        Cl_guess = Cl_guess2
        error_Cl = error_Cl2
        Cl_guess2 = -error_Cl2/slope + Cl_guess
        ! checking for neg values
        if (Cl_guess2 .lt. 0d0) then
            Cl_guess2 = 0d0     
        end if
        ! we can then get the ionic strength for these values of Na and Cl
        ion_guess = 0.5d0 * Na_needed + 0.5d0 * Cl_guess2
        ! and call the elctroneutrality again
        call calc_electroneutrality_error(ion_guess, Na_needed, Cl_guess2)
        error_Cl2 = electro_error

        ! checking for convergence
        if (error_Cl2 .eq. 0d0) then
            break = 1 ! we have convergence!
        end if
        if (abs(Cl_guess/Cl_guess2 -1d0) .lt. 1d-25) then
            break = 2 ! no more change in the value of x!
        end if
        if (abs(error_Cl/error_Cl2 - 1d0) .lt. 1d-25) then
            break = 3 ! no more change in the value of f(x)
        end if
        if (counter .gt. 1000) then
            break = -1 ! no convergence
        end if
        counter = counter +1
    end do
    get_Cl_conc = Cl_guess2

    end function get_Cl_conc

    double precision function get_ion_str_real (ion_guess, Na_guess, Cl_guess)
    ! gets the real ionic strength!

    ! parameters
    implicit none
    double precision :: ion_guess, Na_guess, Cl_guess
    double precision :: ion_guess2, ion_error, ion_error2, slope_ion
    integer ::  counter

    ! we can then get the ionic strength for these values of Na and Cl
    ion_break= 0 ! initialising
    counter = 1 ! initialising!
    ! initial point 1
    call calc_ion_str_error(ion_guess, Na_guess, Cl_guess)
    ion_error = ion_str_error
    ! and the second initial point
    ion_guess2 = 2d0*ion_guess + 30d-3
    call calc_ion_str_error(ion_guess2, Na_guess, Cl_guess)
    ion_error2 = ion_str_error
    ! we then use a derivative method!
    do while (ion_break .eq. 0)
        slope_ion = (ion_error2-ion_error) / (ion_guess2 - ion_guess)
        ion_guess = ion_guess2
        ion_error = ion_error2
        ion_guess2 = -ion_error2/slope_ion + ion_guess
        call calc_ion_str_error(ion_guess2, Na_guess, Cl_guess)
        ion_error2 = ion_str_error

        ! checking for convergance
        if (ion_error2 .eq. 0d0) then
            ion_break = 1 ! perfect convergence!
        end if
        if (abs(ion_guess/ion_guess2 -1d0) .lt. 1d-25) then
            ion_break =2 ! no more changes in the values of x!
        end if
        if (abs(ion_error/ion_error2 -1d0) .lt. 1d-25) then
            ion_break = 3 ! no more changes in the value of f(x)
        end if
        if (counter .gt. 1000) then
            ion_break = -1 ! convergence not achieved!
            exit
        end if
        counter = counter +1
    end do

    get_ion_str_real = ion_guess2

    end function get_ion_str_real

    subroutine calc_electroneutrality_error(ion_guess, Na_guess, Cl_guess)
    ! gets the electroneutrality error!

    ! parameters
    implicit none
    double precision :: ion_guess, Na_guess, Cl_guess

    ! get the species concentration
    call calc_species_conc(ion_guess)

    ! and the electroneutrality
    electro_error = (-AcO-H2PO4-2d0*HPO4-3d0*PO4-succH-2d0*succ+tris-H2cit-2d0*Hcit &
        -3d0*cit-malaH-2d0*mala-arg1+arg2+arg3-his1+his2+his3-cys1+cys2-cys3 &
        -H3EDTA-2d0*H2EDTA-3d0*HEDTA-4d0*EDTA+H-OH+Na_guess-Cl_guess)



    end subroutine calc_electroneutrality_error

    subroutine calc_ion_str_error (ion_str_guess, Na_guess, Cl_guess)
    ! gets the ion_str error!

    ! parameters
    implicit none
    double precision :: ion_str_guess, ion_str_calc, Na_guess, Cl_guess

    ! getting the species values
    call calc_species_conc(ion_str_guess)

    ! and calculating the error!
    ion_str_calc = 0.5d0*(1d0**2d0 * AcO)&
        +  0.5d0 *( 1d0**2d0 * H2PO4) + 0.5d0 *(2d0**2d0 * HPO4) + 0.5d0*(3d0**2d0 * PO4) &
        + 0.5d0*(1d0**2d0 * succH) + 0.5d0 *(2d0**2d0 * succ) &
        + 0.5d0*(1d0**2d0 * tris) &
        + 0.5d0*(1d0**2d0 * H2cit) + 0.5d0*(2d0**2d0 * Hcit) + 0.5d0*(3d0**2d0 * cit) &
        + 0.5d0*(1d0**2d0 * malaH) + 0.5d0*(2d0**2d0 * mala) &
        + 0.5d0*(1d0**2d0 *arg3) + 0.5d0*(1d0**2d0 * his3) + 0.5d0*(1d0**2d0 * cys3) &
        + 0.5d0*(1d0**2d0 * H3EDTA) + 0.5d0*(2d0**2d0 * H2EDTA) + 0.5d0*(3d0**2d0 * HEDTA) + 0.5d0*(4d0**2d0 * EDTA) &
        + 0.5d0*(H + OH + Na_guess + Cl_guess)

    ion_str_error = ion_str_calc - ion_str_guess

    end subroutine calc_ion_str_error

    subroutine calc_species_conc(ion_str_guess)
    ! gets the  concentrations of all interesting species!

    ! parameters
    implicit none
    double precision :: pKa_eff(4), ion_str_guess,  denom
    integer :: counter, loop

    ! juste starting some parameters
    do loop = 1,4
        pKa_eff(loop)=0d0
    end do

    H = 10**-in_pH
    OH = (10**-pKw)/H

    ! acetate
    counter = 1
    pKa_eff(1) = pKa_acet + (2*(z_acet-1)*(((A_DH*sqrt(ion_str_guess))/(1+sqrt(ion_str_guess)))-B_DH(counter)*ion_str_guess))
    denom = 1+((10**-pKa_eff(1))/H)
    AcO = in_molarities(counter)*10**(-pKa_eff(1))/H/denom

    ! phosphate
    counter = counter + 1
    pKa_eff(1) = pKa_phos(1) + (2*(z_phos(1)-1)*(((A_DH*sqrt(ion_str_guess))/(1+sqrt(ion_str_guess)))-B_DH(counter)*ion_str_guess))
    pKa_eff(2) = pKa_phos(2) + (2*(z_phos(2)-1)*(((A_DH*sqrt(ion_str_guess))/(1+sqrt(ion_str_guess)))-B_DH(counter)*ion_str_guess))
    pKa_eff(3) = pKa_phos(3) + (2*(z_phos(3)-1)*(((A_DH*sqrt(ion_str_guess))/(1+sqrt(ion_str_guess)))-B_DH(counter)*ion_str_guess))
    denom = 1+((10**-pKa_eff(1))/H)+((10**-pKa_eff(1))*(10**-pKa_eff(2))/H**2d0)+((10**-pKa_eff(1))*(10**-pKa_eff(2))*(10**-pKa_eff(3))/H**3d0)
    H2PO4 = in_molarities(counter)* 10**(-pKa_eff(1))/H/denom
    HPO4 = H2PO4 * 10**(-pKa_eff(2))/H
    PO4 = HPO4 * 10**(-pKa_eff(3))/H

    ! succinate
    counter = counter + 1
    pKa_eff(1) = pKa_succ(1) + (2*(z_succ(1)-1)*(((A_DH*sqrt(ion_str_guess))/(1+sqrt(ion_str_guess)))-B_DH(counter)*ion_str_guess))
    pKa_eff(2) = pKa_succ(2) + (2*(z_succ(2)-1)*(((A_DH*sqrt(ion_str_guess))/(1+sqrt(ion_str_guess)))-B_DH(counter)*ion_str_guess))
    denom = 1+((10**-pKa_eff(1))/H)+((10**-pKa_eff(1))*(10**-pKa_eff(2))/H**2d0)
    succH = in_molarities(counter)* 10**(-pKa_eff(1))/H/denom
    succ= succH * 10**(-pKa_eff(2))/H


    ! tris
    counter = counter + 1
    pKa_eff(1) = pKa_tris + ((2*z_tris-1)*(((A_DH*sqrt(ion_str_guess))/(1+sqrt(ion_str_guess)))-B_DH(counter)*ion_str_guess))
    denom = 1+((10**-pKa_eff(1))/H)
    tris = in_molarities(counter) - in_molarities(counter)* 10**(-pKa_eff(1))/H/denom

    ! citrate
    counter = counter + 1
    pKa_eff(1) = pKa_citr(1) + (2*(z_citr(1)-1)*(((A_DH*sqrt(ion_str_guess))/(1+sqrt(ion_str_guess)))-B_DH(counter)*ion_str_guess))
    pKa_eff(2) = pKa_citr(2) + (2*(z_citr(2)-1)*(((A_DH*sqrt(ion_str_guess))/(1+sqrt(ion_str_guess)))-B_DH(counter)*ion_str_guess))
    pKa_eff(3) = pKa_citr(3) + (2*(z_citr(3)-1)*(((A_DH*sqrt(ion_str_guess))/(1+sqrt(ion_str_guess)))-B_DH(counter)*ion_str_guess))
    denom = 1+((10**-pKa_eff(1))/H)+((10**-pKa_eff(1))*(10**-pKa_eff(2))/H**2d0)+((10**-pKa_eff(1))*(10**-pKa_eff(2))*(10**-pKa_eff(3))/H**3d0)
    H2cit = in_molarities(counter)* 10**(-pKa_eff(1))/H/denom
    Hcit = H2cit * 10**(-pKa_eff(2))/H
    cit = Hcit * 10**(-pKa_eff(3))/H

    ! malate
    counter = counter + 1
    pKa_eff(1) = pKa_mala(1) + (2*(z_mala(1)-1)*(((A_DH*sqrt(ion_str_guess))/(1+sqrt(ion_str_guess)))-B_DH(counter)*ion_str_guess))
    pKa_eff(2) = pKa_mala(2) + (2*(z_mala(2)-1)*(((A_DH*sqrt(ion_str_guess))/(1+sqrt(ion_str_guess)))-B_DH(counter)*ion_str_guess))
    denom = 1+((10**-pKa_eff(1))/H)+((10**-pKa_eff(1))*(10**-pKa_eff(2))/H**2d0)
    malaH = in_molarities(counter)* 10**(-pKa_eff(1))/H/denom
    mala= malaH * 10**(-pKa_eff(2))/H

    ! argineine
    ! arg1
    counter = counter + 1
    pKa_eff(1) = pKa_argi(1) + (2*(z_argi(1)-1)*(((A_DH*sqrt(ion_str_guess))/(1+sqrt(ion_str_guess)))-B_DH(counter)*ion_str_guess))
    denom = 1+((10**-pKa_eff(1))/H)
    arg1 = in_molarities(counter)*10**(-pKa_eff(1))/H/denom
    ! arg2
    pKa_eff(2) = pKa_argi(2) + ((2*z_argi(2)-1)*(((A_DH*sqrt(ion_str_guess))/(1+sqrt(ion_str_guess)))-B_DH(counter)*ion_str_guess))
    denom = 1+((10**-pKa_eff(2))/H)
    arg2 = in_molarities(counter) - in_molarities(counter)* 10**(-pKa_eff(2))/H/denom
    ! arg3
    pKa_eff(3) = pKa_argi(3) + ((2*z_argi(3)-1)*(((A_DH*sqrt(ion_str_guess))/(1+sqrt(ion_str_guess)))-B_DH(counter)*ion_str_guess))
    denom = 1+((10**-pKa_eff(3))/H)
    arg3 = in_molarities(counter) - in_molarities(counter)* 10**(-pKa_eff(3))/H/denom

    ! histidine
    ! his1
    counter = counter + 1
    pKa_eff(1) = pKa_hist(1) + (2*(z_hist(1)-1)*(((A_DH*sqrt(ion_str_guess))/(1+sqrt(ion_str_guess)))-B_DH(counter)*ion_str_guess))
    denom = 1+((10**-pKa_eff(1))/H)
    his1 = in_molarities(counter)*10**(-pKa_eff(1))/H/denom
    ! his2
    pKa_eff(2) = pKa_hist(2) + ((2*z_hist(2)-1)*(((A_DH*sqrt(ion_str_guess))/(1+sqrt(ion_str_guess)))-B_DH(counter)*ion_str_guess))
    denom = 1+((10**-pKa_eff(2))/H)
    his2 = in_molarities(counter) - in_molarities(counter)* 10**(-pKa_eff(2))/H/denom
    ! his3
    pKa_eff(3) = pKa_hist(3) + ((2*z_hist(3)-1)*(((A_DH*sqrt(ion_str_guess))/(1+sqrt(ion_str_guess)))-B_DH(counter)*ion_str_guess))
    denom = 1+((10**-pKa_eff(3))/H)
    his3 = in_molarities(counter) - in_molarities(counter)* 10**(-pKa_eff(3))/H/denom

    ! cysteine
    ! cys1
    counter = counter + 1
    pKa_eff(1) = pKa_cyst(1) + (2*(z_cyst(1)-1)*(((A_DH*sqrt(ion_str_guess))/(1+sqrt(ion_str_guess)))-B_DH(counter)*ion_str_guess))
    denom = 1+((10**-pKa_eff(1))/H)
    cys1 = in_molarities(counter)*10**(-pKa_eff(1))/H/denom
    ! cys2
    pKa_eff(2) = pKa_cyst(2) + ((2*z_cyst(2)-1)*(((A_DH*sqrt(ion_str_guess))/(1+sqrt(ion_str_guess)))-B_DH(counter)*ion_str_guess))
    denom = 1+((10**-pKa_eff(2))/H)
    cys2 = in_molarities(counter) - in_molarities(counter)* 10**(-pKa_eff(2))/H/denom
    ! cys3
    pKa_eff(3) = pKa_cyst(3) + (2*(z_cyst(3)-1)*(((A_DH*sqrt(ion_str_guess))/(1+sqrt(ion_str_guess)))-B_DH(counter)*ion_str_guess))
    denom = 1+((10**-pKa_eff(3))/H)
    cys3 = in_molarities(counter)*10**(-pKa_eff(3))/H/denom

    ! EDTA
    counter = counter + 1
    pKa_eff(1) = pKa_EDTA(1) + (2*(z_EDTA(1)-1)*(((A_DH*sqrt(ion_str_guess))/(1+sqrt(ion_str_guess)))-B_DH(counter)*ion_str_guess))
    pKa_eff(2) = pKa_EDTA(2) + (2*(z_EDTA(2)-1)*(((A_DH*sqrt(ion_str_guess))/(1+sqrt(ion_str_guess)))-B_DH(counter)*ion_str_guess))
    pKa_eff(3) = pKa_EDTA(3) + (2*(z_EDTA(3)-1)*(((A_DH*sqrt(ion_str_guess))/(1+sqrt(ion_str_guess)))-B_DH(counter)*ion_str_guess))
    pKa_eff(4) = pKa_EDTA(4) + (2*(z_EDTA(4)-1)*(((A_DH*sqrt(ion_str_guess))/(1+sqrt(ion_str_guess)))-B_DH(counter)*ion_str_guess))
    denom = 1+((10**-pKa_eff(1))/H)+((10**-pKa_eff(1))*(10**-pKa_eff(2))/H**2d0)+((10**-pKa_eff(1))*(10**-pKa_eff(2))*(10**-pKa_eff(3))/H**3d0)+((10**-pKa_eff(1))*(10**-pKa_eff(2))*(10**-pKa_eff(3))*(10**-pKa_eff(4))/H**4d0)
    H3EDTA = in_molarities(counter)* 10**(-pKa_eff(1))/H/denom
    H2EDTA = H3EDTA * 10**(-pKa_eff(2))/H
    HEDTA = H2EDTA * 10**(-pKa_eff(3))/H
    EDTA = HEDTA * 10**(-pKa_eff(4))/H

    end subroutine calc_species_conc

    end module calc_grams

