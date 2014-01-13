C
C File:   testflux.for
C Author: ash
C
C Created on 13 Январь 2014 г., 14:20
C
      LOGICAL :: Flux_init,Flux_open_file,Flux_set_sbit
      LOGICAL :: Flux_read_hdr,Flux_print_table
      LOGICAL :: Flux_close_file
      logical res,tmp
      character(*),parameter:: FileName="test.flux            "
      PRINT*, 'Hello World'
      res=Flux_init()
      res=Flux_set_sbit(1,.TRUE.)
      res=Flux_set_sbit(3,.TRUE.)
      res=Flux_open_file(FileName)
      do while(res)
          res=Flux_read_hdr()
          if(res)tmp=Flux_print_table()
      end do
      res=Flux_close_file()
      END

