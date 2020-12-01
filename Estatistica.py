import numpy as np
from scipy import math
from scipy.stats import norm, t, f, chi2, f_oneway
class Stat:
    def __init__(self, A=None, B=None, alfa=0.05):
        self.A = A
        self.B = B
        self.alfa = alfa
        self.calcStat()

    def calcStat(self):
        if self.A:
            self.Na = len(self.A)
            self.Xa = np.mean(self.A)
            self.Va = np.var(self.A)
        if self.B:
            self.Nb = len(self.B)       
            self.Xb = np.mean(self.B)
            self.Vb = np.var(self.B)
        return None

    def setData(self, A=None, B=None):
        if A: self.A = A
        if B: self.B = B
        self.calcStat()
    
    def setData2(self, Xa = None, Xb = None, Va = None, Vb = None, Na = None, Nb=None, Pla = None, Plb = None):
        if Xa: self.Xa = Xa
        if Xb: self.Xb = Xb
        if Va: self.Va = Va
        if Vb: self.Vb = Vb
        if Na: self.Na = Na
        if Nb: self.Nb = Nb
        if Pla: self.pla = Pla
        if Plb: self.plb = Plb
        self.calcStat()

    def setConf(self, alfa):
        self.alfa = alfa

    def getData(self, data='AB'):
        if data == 'A' or data =='a': return self.A
        if data == 'B' or data =='b': return self.B
        if data == 'AB' or data =='ab': return [self.A, self.B]

    def getConf(self):
        return self.alfa

    def IC(self, mode = "diffMed", student = False):
        '''Calcula o intervalo de confinaça para um determinado conjunto de 
        dados'''
        if mode == 'diffMed':
            if student == False: Zcrit = norm.ppf(1-self.alfa/2)
            else: Zcrit = t.ppf(1-self.alfa/2)
            
            e0 = Zcrit * np.sqrt(self.Va / self.Na + self.Vb / self.Nb)
            dX = self.Xa - self.Xb
            print("IC = ", int(1000*(dX))/1000, " ± ", int(1000*(e0))/1000)
            return (int(1000*(dX - e0))/1000, int(1000*(dX + e0))/1000)

        if mode == 'diffProp':
            if self.pla>1: self.pla /= self.Na
            if self.plb>1: self.plb /= self.Nb

            if student == False: Zcrit = norm.ppf(1-self.alfa/2)
            else: Zcrit = t.ppf(1-self.alfa/2)

            e0 = Zcrit * np.sqrt(self.pla * (1-self.pla) / self.Na + self.plb * (1-self.plb) / self.Nb)
            dX = self.pla - self.plb

            print("IC = ", int(1000*(dX))/1000, " ± ", int(1000*(e0))/1000)
            return (int(1000*(dX-e0))/1000, int(1000*(dX+e0))/1000)

        else:
            print("Erro de entrada")

    def TH(self, mode = "med", cond = None, student = False, mesmaVariancia = False):
        '''Realiza um teste de hipóteses dadas as condições de entrada
            Mode: escolha do tipo de analise que esta sendo feita
            Cond: Condição inicial da comparação (por exemplo: mi=25)
            '''

        if mode == "med":
            Zcalc = (self.Xa - cond) / np.sqrt(self.Va / self.Na)
            print("Zcalc = ", int(1000 * Zcalc)/1000)
            
            if student == False:
                Zcritmin = norm.ppf(self.alfa/2)
                Zcritmax = norm.ppf(1- self.alfa/2)
            else:
                Zcritmin = t.ppf(self.alfa/2, self.Na - 1)
                Zcritmax = t.ppf(1- self.alfa/2, self.Na - 1)
            
            print("Zcritmin = ", int(1000 * Zcritmin)/1000)
            print("Zcritmax = ", int(1000 * Zcritmax)/1000)

            if Zcritmin <= Zcalc <= Zcritmax:
                print("Aceito H0")
                return True
            else:
                print("Rejeito H0")
                return False

        elif mode == 'diffMedEmp' or mode == "maiorMedEmp" or mode == "menorMedEmp":
            # Modo para diferença de dados emparelhados
            if mesmaVariancia == False: 
                dif = [self.A[i] - self.B[i] for i in range(self.Na)]
                sd2 = np.var(dif, ddof=1)
                T_calc_emparelhado = np.mean(dif) * np.sqrt(len(dif) / sd2)
                grausLib = self.Na
            else:
                sd2 = self.Va / self.Na + self.Vb / self.Nb
                T_calc_emparelhado = (self.Xa - self.Xb) / np.sqrt(sd2)
                grausLib = self.Na + self.Nb

            print("T_calc_emparelhado = ", int(1000 * T_calc_emparelhado)/1000)

            if student == False: 
                return("Ainda nao implementado. Por favor, utilize a T-Student")

            else: 
                if mode == "diffMedEmp": 
                    Tcritmin = t.ppf(self.alfa/2, grausLib-1)
                    Tcritmax = t.ppf(1-self.alfa/2, grausLib-1)
                elif mode == "maiorMedEmp": Tcrit = t.ppf(1-self.alfa, grausLib-1)
                else: Tcrit = t.ppf(self.alfa, grausLib-1)

                if mode == "diffMedEmp":
                    if Tcritmin < T_calc_emparelhado < Tcritmax:
                        print("Tcritmin = ", int(1000 * Tcritmin)/1000)
                        print("Tcritmax = ", int(1000 * Tcritmax)/1000)
                        print("Aceito H0, as médias são iguais")
                        print("Lembre que deve ser respeitada a ordem de que o dado A é maior que o B. A troca de sinal irá distorcer o resultado")
                        return True
                    else:
                        print("Tcritmin = ", int(1000 * Tcritmin)/1000)
                        print("Tcritmax = ", int(1000 * Tcritmax)/1000)
                        print("Rejeito H0, as médias são diferentes")
                        return False
                
                elif mode == "maiorMedEmp":
                    if T_calc_emparelhado > Tcrit:
                        print("Tcrit = ", int(1000 * Tcrit)/1000)
                        print("Rejeito H0, a média A é maior que a média B")
                        print("Lembre que deve ser respeitada a ordem de que o dado A é menor que o B. A troca de sinal irá distorcer o resultado")
                        return False
                    else:
                        print("Tcrit = ", int(1000 * Tcrit)/1000)
                        print("Aceito H0, as médias são iguais")
                        print("Lembre que deve ser respeitada a ordem de que o dado A é menor que o B. A troca de sinal irá distorcer o resultado")
                        return True

                else:
                    if T_calc_emparelhado < Tcrit:
                        print("Tcrit = ", int(1000 * Tcrit)/1000)
                        print("Rejeito H0, a média A é menor que a média B")
                        return False
                    else:
                        print("Tcrit = ", int(1000 * Tcrit)/1000)
                        print("Aceito H0, as médias são iguais")
                        return True
                   

        elif mode == "menorVar":
            Quicalc = (self.Na - 1) * self.Va / cond
            Quicrit = chi2.ppf(self.alfa, self.Na-1)
            print("Quicalc = ", int(1000*(Quicalc))/1000)
            print("Quicritmin = ", int(1000*(Quicrit))/1000)
            if Quicalc < Quicrit:
                print("Rejeito H0")
                return False
            else:
                print("Aceito H0")
                return True

        elif mode == "maiorVar":
            Quicalc = (self.Na - 1) * self.Va / cond
            Quicrit = chi2.ppf(1-self.alfa, self.Na-1)
            print("Quicalc = ", int(1000*(Quicalc))/1000)
            print("Quicritmax = ", int(1000*(Quicrit))/1000)
            if Quicalc < Quicrit:
                print("Aceito H0")
                return True
            else:
                print("Rejeito H0")
                return False

        elif mode == "diffVar":
            Fcalc = self.Va/ self.Vb
            Fcritmin = f.ppf(self.alfa/2, self.Na-1, self.Nb-1)
            Fcritmax = f.ppf(1-self.alfa/2, self.Na-1, self.Nb-1)
            print("Fcalc = ", int(1000*(Fcalc))/1000)
            print("Fcritmin = ", int(1000*(Fcritmin))/1000)
            print("Fcritmax = ", int(1000*(Fcritmax))/1000)

            if Fcritmin <= Fcalc <= Fcritmax:
                print("Aceito H0")
                return True
            else:
                print("Rejeito H0")
                return False
    def pValor(self, mode = "media", cond = None):
        '''Calcula o p-valor dadas as condições de entrada
            Mode: escolha do tipo de analise que esta sendo feita
            Cond: Condição inicial da comparação (por exemplo: mi=25)
            '''

        if mode == "media":
            Zcalc = (self.Xa - cond) / np.sqrt(self.Va / self.Na)
            temp = norm.cdf(Zcalc)

            print("Alfa = ", int(1000*(100 * (1 - temp)))/1000, "%")
            return 1 - temp

        if mode == "alfaMedia":
            '''A diferença deste método para o "media" é que o "cond" é um critério. Ele serve para calcular o erro alfa quando
            se tem um conjunto de dados de origem desconhecida e adota-se um valor (arbitrátrio) para decidir se ele pertence à
            um certo conjunto de dados A ou certo conjunto de dados B. Vide exercicio 9 da P2 de estatística de 2018'''

            Zcalc = (cond - self.Xa) / np.sqrt(self.Va / self.Na)
            temp = norm.cdf(Zcalc)

            print("Alfa = {:.3f}".format(100 * (1 - temp)), "%")
            return 1 - temp

        elif mode == "betaMedia":
            Zcalc = (cond - self.Xb) / np.sqrt(self.Vb / self.Nb)
            temp = norm.cdf(Zcalc)

            print("Beta = {:.3f}".format(100 * temp), "%")
            return temp

        elif mode == "maiorProp":
            # Apenas para tornar o codigo flexivel à entradas com o valor
            # inteiro ao invés de inserir a proporção
            if self.pla>1 and self.pla: self.pla /= self.Na
            if cond > 1: cond /= self.Na

            Zcalc = (self.pla - cond) / np.sqrt(cond * (1 - cond) / self.Na)
            temp = norm.cdf(Zcalc)

            print("Alfa = ", int(1000*(1000 * (1 - temp)))/1000, "%")
            return 1-temp

        else:
            print("Erro de entrada")


    def calculateTA(self, A, B):
        temp = 0
        N = 0
        for i in range(len(A)):
            temp += B[i] * A[i]
            N += B[i]

        return temp / N, N

    def TA(self, A, B, mode = "poisson"):
        '''Realiza um teste de aderência dadas as condições de entrada
        Para o teste, A deve conter o número de trocas, defeitos, etc, e 
        B deve conter quantas vezez essa troca ou defeito ocorreu. Note que 
        é fundamental que haja compatibilidade perfeita entre as entradas'''

        if len(A) != len(B): return print("Entradas incompatíveis")
        X, N = self.calculateTA(A, B)

        if mode == "poisson":
            Oi = B
            Chi2calc = 0
            n = 0
            Oiacc = 0
            Eiacc = 0

            # Generate the Expected number by the poisson model
            for i in range(len(A)):
                Eiacc += np.exp(-X) * X**A[i] / math.factorial(A[i]) * N
                Oiacc += Oi[i]
                
                if i != (len(A)-1) and np.exp(-X) * X**A[i+1] / math.factorial(A[i+1]) * N >= 5:
                    Chi2calc += (Oiacc - Eiacc)**2 / Eiacc
                    n += 1
                    Oiacc = 0
                    Eiacc = 0

            Chi2calc += (Oiacc - Eiacc)**2 / Eiacc
            n += 1  
                
                
            Chi2crit = chi2.ppf(1-self.alfa, n - 2)
        
            print("Chi2calc = ", int(1000 * Chi2calc) / 1000)
            print("Chi2crit = ", int(1000 * Chi2crit) / 1000)
            if Chi2calc < Chi2crit:
                print("Aceito H0")
                return True
            else:
                print("Rejeito H0")
                return False

    def AV(self, SQ1 = None, n1 = None, SQ2 = None, n2 = None,  SQR = None, SQT = None, QM1 = None, QM2 = None, QMR = None, QMT = None, mode = 'multiMedia'):
        if mode == 'multiMedia':
            gl1 = n1 - 1
            glr = n1 * (n2-1)
            glt = n1 * n2 - 1

            if QM1:
                SQ1 = QM1 * gl1
            if QMR:
                SQR = QMR * glr
            if QMT:
                SQT = QMT * glt
            if not SQ1:
                SQ1 = SQT - SQR
            if not SQR:
                SQR = SQT - SQ1
            if not SQT:
                SQT = SQ1 + SQR

            QM1 = SQ1 / gl1
            QMR = SQR / glr
            QMT = SQT / glt

            F1calc = QM1 / QMR
            F1crit = f.ppf(1 - self.alfa, gl1, glr)
            pVal = f.cdf(F1calc, gl1, glr)


            print('| Fonte        SQ        gl          QM        F    Fcrit    pVal  |')
            print('| Dado 1    {:.2f}     {:.0f}       {:.3f}   {:.3f}  {:.3f}  {:.5f} |'.format(SQ1, gl1, QM1, F1calc, F1crit, pVal))
            print('| Residual  {:.2f}    {:.0f}      {:.3f} |'.format(SQR, glr, QMR))
            print('| Total     {:.2f}    {:.0f}      {:.3f} |'.format(SQT, glt, QMT))
            print()

            if F1calc < F1crit:
                print("Aceito H0, as médias são iguais")
                return True
            else:
                print("Rejeito H0, as médias são diferentes")
                return False


        elif mode == 'doisFat':
            gl1 = n1 - 1
            gl2 = n2 - 1
            glr = gl1 * gl2
            glt = n1 * n2 - 1

            if QM1:
                SQ1 = QM1 * gl1
            if QM2:
                SQ2 = QM2 * gl2
            if QMR:
                SQR = QMR * glr
            if QMT:
                SQT = QMT * glt
            if not SQ1:
                SQ1 = SQT - SQR - SQ2
            elif not SQ2:
                SQ2 = SQT - SQR - SQ1
            elif not SQR:
                SQR = SQT - SQ1 - SQ2
            elif not SQT:
                SQT = SQ1 + SQ2 + SQR

            QM1 = SQ1 / gl1
            QM2 = SQ2 / gl2
            QMR = SQR / glr
            QMT = SQT / glt

            F1calc = QM1 / QMR
            F2calc = QM2 / QMR

            F1crit = f.ppf(1 - self.alfa, gl1, glr)
            F2crit = f.ppf(1 - self.alfa, gl2, glr)

            print('| Fonte       SQ        gl        QM      F     Fcrit |')
            print('| Dado 1    {:.2f}       {:.0f}       {:.3f}  {:.3f}  {:.3f} |'.format(SQ1, gl1, QM1, F1calc, F1crit))
            print('| Dado 2    {:.2f}       {:.0f}       {:.3f}   {:.3f}   {:.3f} |'.format(SQ2, gl2, QM2, F2calc, F2crit))
            print('| Residual  {:.2f}       {:.0f}      {:.3f} |'.format(SQR, glr, QMR))
            print('| Total     {:.2f}      {:.0f}      {:.3f} |'.format(SQT, glt, QMT))
            print()

            if F1calc < F1crit:
                print("Aceito H0.1, as médias de 1 são iguais")
                R1 = True
            else:
                print("Rejeito H0.1, as médias de 1 são diferentes")
                R1 = False
            
            if F2calc < F2crit:
                print("Aceito H0.2, as médias de 2 são iguais")
                return (R1, True)
            else:
                print("Rejeito H0.2, as médias de 2 são diferentes")
                return (R1, False)