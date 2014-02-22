!**********************************************************************!
ENTRY Flux_Get_Emin(iNuAnu,iFlavor)
!----------------------------------------------------------------------!
    Flux_Get_Emin=E_min(iNuAnu,iFlavor)
    return

!**********************************************************************!
ENTRY Flux_Get_Emax(iNuAnu,iFlavor)
!----------------------------------------------------------------------!
    Flux_Get_Emax=E_max(iNuAnu,iFlavor)
    return

!**********************************************************************!
ENTRY Flux_Get_Zmin(iNuAnu,iFlavor)
!----------------------------------------------------------------------!
    Flux_Get_Zmin=Z_min(iNuAnu,iFlavor)
    return

!**********************************************************************!
ENTRY Flux_Get_Zmax(iNuAnu,iFlavor)
!----------------------------------------------------------------------!
    Flux_Get_Zmax=Z_max(iNuAnu,iFlavor)
    return

!**********************************************************************!
ENTRY Flux_Get_LastNu(iNuAnu,iFlavor)                                  !do we really need this?
!----------------------------------------------------------------------!
    write(*,*) NuAnu,Flavor
    iNuAnu=NuAnu; iFlavor=Flavor
    return

!**********************************************************************!
ENTRY Flux_Open_File(file_name)
!----------------------------------------------------------------------!
!space removing from file_name-----------------------------------------!
    filename=trim(file_name)                                           !does it really work? it seems no
!----------------------------------------------------------------------!by the way we can use C-string for this!
    write(*,*) 'Opening of flux file: ',filename
    Flux_Open_File=.false.
    open(flxf,status='old',file=filename,err=98)
    Flux_Open_File=.true.
    return

!**********************************************************************!
ENTRY Flux_Read_Hdr()
!----------------------------------------------------------------------!
    LogE=.true.                                                        !Default
!table header reading with printing reports----------------------------!
    write(*,*) 'Table header:'
    write(*,'(A29)') hf
    do
        read(flxf,*,iostat=ioer) fline
!testing: end of file--------------------------------------------------!
        if(ioer.ne.0)then
            Flux_Read_Hdr=.false.
            write(*,*) 'EOF reached!'
            return
        endif
!testing: not a header-------------------------------------------------!
        if(fline(1:1).ne.'#')exit
!table header line reading---------------------------------------------!
        if(fline(2:7)=='Ntype=')then
            read(fline,'(7X,A2)') NuType
            write(*,'(A1,1X,"Neutrino type:",9X,A2,1X,A1)') vf,NuType,vf
            selectcase(NuType(1:1))
                case('n'); NuAnu=1
                case('a'); NuAnu=2
                case default
                    Flux_Read_Hdr=.false.
                    write(*,*) 'Flux_Read_Hdr ERROR: Unknown NuAnu: ',NuType(1:1),'!'
                    return
            end select
            selectcase(NuType(2:2))
                case('e'); Flavor=1
                case('m'); Flavor=2
                case('t'); Flavor=3
                case default
                    Flux_Read_Hdr=.false.
                    write(*,*) 'Flux_Read_Hdr ERROR: Unknown Flavor: ',NuType(2:2),'!'
                    return
            endselect
        elseif(fline(2:4)=='NE=')then
            read(fline,'(4X,I5)') NE(NuAnu,Flavor)
            write(*,'(A1,1X,"Energy point number:",1X,I4,1X,A1)') vf,NE(NuAnu,Flavor),vf
        elseif(fline(2:6)=='NZ=')then
            read(fline,'(6X,I5)') NZ(NuAnu,Flavor)
            write(*,'(A1,1X,"Cosine point number:",1X,I4,1X,A1)') vf,NZ(NuAnu,Flavor),vf
        elseif(fline(2:5)=='RegE')then
            LogE=.false.
            write(*,'(A1,1X,"Energy scale:",5X,A7,1X,A1)') vf,'regular',vf
        elseif(fline(2:5)=='LogE')then
            LogE=.true.
            write(*,'(A1,1X,"Energy scale:",1X,A11,1X,A1)') vf,'logarithmic',vf
!----------------------------------------------------------------------!
        else
            write(*,*) 'Unknown line: "',fline,'"!'
        endif
    enddo
    write(*,'(A29)') hf
!----------------------------------------------------------------------!
    HasTable(NuAnu,Flavor)=.true.
    LogTable(NuAnu,Flavor)=LogE
    Flux_Read_Hdr=.true.
    return

!**********************************************************************!
ENTRY Flux_Read_Table()
!----------------------------------------------------------------------!
!flux table reading----------------------------------------------------!
    write(*,*) 'Flux table reading...'
    do n_E=1,NE(NuAnu,Flavor)
        read(flxf,*,iostat=ioer) E(NuAnu,Flavor,n_E),dF(NuAnu,Flavor,n_E)
!testing: unexpected end of table--------------------------------------!
        if(ioer.ne.0)then
            Flux_Read_Table=.false.
            stop 'Flux_Read_Table ERROR: unexpected end of flux table!'
        endif
!----------------------------------------------------------------------!
    enddo
!energy limit getting--------------------------------------------------!
    E_min(NuAnu,Flavor)=E(NuAnu,Flavor,1)
    E_max(NuAnu,Flavor)=E(NuAnu,Flavor,NE(NuAnu,Flavor))
    if(E_max(NuAnu,Flavor)<E_min(NuAnu,Flavor))then
        E_min(NuAnu,Flavor)=E_max(NuAnu,Flavor)
        E_max(NuAnu,Flavor)=E(NuAnu,Flavor,1)
    endif
!----------------------------------------------------------------------!
    Flux_Read_Table=.true.
    return

!**********************************************************************!
ENTRY Flux_Close_File()
!----------------------------------------------------------------------!
    close(flxf)
    write(*,*) 'Flux file is closed'
    return

!**********************************************************************!
ENTRY Flux_Has_Table(iNuAnu,iFlavor)
!----------------------------------------------------------------------!
    Flux_Has_Table=HasTable(iNuAnu,iFlavor)
    return

!**********************************************************************!
ENTRY Flux_Has_Spline(iNuAnu,iFlavor)
!----------------------------------------------------------------------!
    Flux_Has_Spline=HasSpline(iNuAnu,iFlavor)
    return

!**********************************************************************!
ENTRY Flux_Print_Table(iNuAnu,iFlavor)
!----------------------------------------------------------------------!
!table property printing-----------------------------------------------!
    write(*,'(A29)') hf
    write(*,'(A1,1X,"Neutrino type set to:",2(1X,I1),1X,A1)') vf,iNuAnu,iFlavor,vf
    write(*,'(A1,1X,"Table size:",5X,I4,"x",I4,1X,A1)') vf,NE(iNuAnu,iFlavor),NZ(iNuAnu,iFlavor),vf
    if(LogTable(iNuAnu,iFlavor))then 
        write(*,'(A1,1X,"Energy scale:",1X,A11,1X,A1)') vf,'logarithmic',vf
    else
        write(*,'(A1,1X,"Energy scale:",5X,A7,1X,A1)') vf,'regular',vf
    endif
    write(*,'(A29)') hf
!flux table printing---------------------------------------------------!
    do n_E=1,NE(iNuAnu,iFlavor)
        write(*,'(A1,2(1X,E12.5),1X,A1)') vf,E(iNuAnu,iFlavor,n_E),dF(iNuAnu,iFlavor,n_E),vf
    enddo
    write(*,'(A29)') hf
!----------------------------------------------------------------------!
    return
