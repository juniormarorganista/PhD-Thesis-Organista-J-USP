cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                subroutines subroutines                         c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine exact_solutions(t)
      
      implicit none
      include 'par.nn'
      include 'comm.vare'
      real*8 t
      
      SELECT CASE (ty_mms)
         !!! Case 1: 
         CASE (1)
            !write(*,*) 'Cases select: ', t
            call exact_solutions_case_01(t)
         !!! Case 2: 
         CASE (2)
            call exact_solutions_case_02(t)
         !!! Case 3: 
         CASE (3)
            call exact_solutions_case_03(t)
         !!! Case 4: 
         CASE (4)
            call exact_solutions_case_04(t)
         !!! Case 5: 
         CASE (5)
            call exact_solutions_case_05(t)
         !!! Case 6: 
         CASE (6)
            call exact_solutions_case_06(t)
         !!! Case 7: 
         CASE (7)
            call exact_solutions_case_07(t)
         !!! Case 8: 
         CASE (8)
            call exact_solutions_case_08(t)
         !!! Case 9: 
         CASE (9)
            call exact_solutions_case_09(t)
         !!! Case 10: 
         CASE (10)
            call exact_solutions_case_10(t)
         !!! Case 11: 
         CASE (11)
            call exact_solutions_case_11(t)
         !!! Case 12: 
         CASE (12)
            call exact_solutions_case_12(t)
         !!! Case 13: 
         CASE (13)
            call exact_solutions_case_13(t)
         !!! Case 14: 
         CASE (14)
            call exact_solutions_case_14(t)
      END SELECT

      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine source_term(t)

      implicit none
      include 'par.nn'
      include 'comm.var'
      real*8 t

      SELECT CASE (ty_mms)
         !!! Case 1: 
         CASE (1)
            call source_term_case_01(t)
         !!! Case 2: 
         CASE (2)
            call source_term_case_02(t)
         !!! Case 3: 
         CASE (3)
            call source_term_case_03(t)
         !!! Case 4: 
         CASE (4)
            call source_term_case_04(t)
         !!! Case 5: 
         CASE (5)
            call source_term_case_05(t)
         !!! Case 6: 
         CASE (6)
            call source_term_case_06(t)
         !!! Case 7: 
         CASE (7)
            call source_term_case_07(t)
         !!! Case 8: 
         CASE (8)
            call source_term_case_08(t)
         !!! Case 9: 
         CASE (9)
            call source_term_case_09(t)
         !!! Case 10: 
         CASE (10)
            call source_term_case_10(t)
         !!! Case 11: 
         CASE (11)
            call source_term_case_11(t)
         !!! Case 12: 
         CASE (12)
            call source_term_case_12(t)
         !!! Case 13: 
         CASE (13)
            call source_term_case_13(t)
         !!! Case 14: 
         CASE (14)
            call source_term_case_14(t)
      END SELECT

      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine boundaries_condictions(t)

      implicit none
      include 'par.nn'
      include 'comm.var'
      include 'comm.vare'
      real*8 t
      
      SELECT CASE (ty_mms)
         !!! Case 1: 
         CASE (1)
           call boundaries_condictions_case_01(t)
         !!! Case 2: 
         CASE (2)
           call boundaries_condictions_case_02(t)
         !!! Case 3: 
         CASE (3)
           call boundaries_condictions_case_03(t)
         !!! Case 4: 
         CASE (4)
           call boundaries_condictions_case_04(t)
         !!! Case 5: 
         CASE (5)
           call boundaries_condictions_case_05(t)
         !!! Case 6: 
         CASE (6)
           call boundaries_condictions_case_06(t)
         !!! Case 7: 
         CASE (7)
           call boundaries_condictions_case_07(t)
         !!! Case 8: 
         CASE (8)
           call boundaries_condictions_case_08(t)
         !!! Case 9: 
         CASE (9)
           call boundaries_condictions_case_09(t)
         !!! Case 10: 
         CASE (10)
           call boundaries_condictions_case_10(t)
         !!! Case 11: 
         CASE (11)
           call boundaries_condictions_case_11(t)
         !!! Case 12: 
         CASE (12)
           call boundaries_condictions_case_12(t)
         !!! Case 13: 
         CASE (13)
           call boundaries_condictions_case_13(t)
         !!! Case 14: 
         CASE (14)
           call boundaries_condictions_case_14(t)
      END SELECT

      return
      end
