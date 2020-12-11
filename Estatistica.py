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
        elif mode == "diffMed":
            Zcalc = (self.Xa - self.Xb) / np.sqrt(self.Va / self.Na + self.Vb / self.Nb)
            print("Zcalc = {:.3f}".format(Zcalc))
            
            if student == False:
                Zcritmin = norm.ppf(self.alfa/2)
                Zcritmax = norm.ppf(1- self.alfa/2)
            else:
                Zcritmin = t.ppf(self.alfa/2, self.Na + self.Nb - 1)
                Zcritmax = t.ppf(1- self.alfa/2, self.Na + self.Nb - 1)
            
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
            pValor = chi2.cdf(Quicalc, self.Na -1)
            print("Quicalc = ", int(1000*(Quicalc))/1000)
            print("Quicritmax = ", int(1000*(Quicrit))/1000)
            print("P-Valor = {:.3f}".format(1-pValor))
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
        elif mode == "diffProp":
            if self.pla>1: self.pla /= self.Na
            if self.plb>1: self.plb /= self.Nb

            if student == False: Zcrit = norm.ppf(1-self.alfa/2)
            else: Zcrit = t.ppf(1-self.alfa/2) 

            pl = (self.Na * self.pla + self.Nb * self.plb) / (self.Na + self.Nb)
            Var = pl * (1-pl) * (1/ self.Na + 1 / self.Nb)
            Zcalc = (self.pla - self.plb)/np.sqrt(Var)
            pValor = 1 - norm.cdf(Zcalc)
            
            print("Zcalc = {:.3f}".format(Zcalc))
            print("Zcrit = {:.3f}".format(Zcrit))
            print("P-Valor = {:.5f}".format(pValor))
            if -Zcrit <= Zcalc <= Zcrit:
                print("Aceito H0, as proporções são iguais")
                return True
            else:
                print("Rejeito H0, as proporções não são iguais")
                return False
    def pValor(self, mode = "media", cond = None, student = False):
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

            if student == False: temp = norm.cdf(Zcalc)
            else: temp = t.cdf(Zcalc, self.Na -1)

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

            print("Alfa = {:.3f} %".format(100*(1-temp)))
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

    def TA(self, A= None, B= None, lamb = None, mode = "poisson"):
        '''Realiza um teste de aderência dadas as condições de entrada
        Para o teste, A deve conter o número de trocas, defeitos, etc, e 
        B deve conter quantas vezez essa troca ou defeito ocorreu. Note que 
        é fundamental que haja compatibilidade perfeita entre as entradas'''
        if self.A: A = self.A
        if self.B: B = self.B
        if len(A) != len(B): return print("Entradas incompatíveis")
        X, N = self.calculateTA(A, B)

        if mode == "poisson":
            Oi = B
            Chi2calc = 0
            n = 0
            Oiacc = 0
            Eiacc = 0

            # Generate the Expected number by the poisson model
            if not lamb:
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

            else:
                for i in range(len(A)):
                    Eiacc += np.exp(-lamb) * lamb**A[i] / math.factorial(A[i]) * N
                    Oiacc += Oi[i]
                    
                    if i != (len(A)-1) and np.exp(-X) * X**A[i+1] / math.factorial(A[i+1]) * N >= 5:
                        Chi2calc += (Oiacc - Eiacc)**2 / Eiacc
                        n += 1
                        Oiacc = 0
                        Eiacc = 0

                Chi2calc += (Oiacc - Eiacc)**2 / Eiacc
                n += 1
            print("Graus de Liberdade = {:.3f}".format(n-2))
                
                
            Chi2crit = chi2.ppf(1-self.alfa, n - 2)
        
            print("Chi2calc = {:.3f}".format(Chi2calc))
            print("Chi2crit = {:.3f}".format(Chi2crit))
            if Chi2calc < Chi2crit:
                print("Aceito H0")
                return True
            else:
                print("Rejeito H0")
                return False

    def AV(self, multiDados = None, SQ1 = None, n1 = None, SQ2 = None, n2 = None,  SQI = None, repeticao = None, SQR = None, SQT = None, QM1 = None, QM2 = None, QMI = None, QMR = None, QMT = None, mode = 'multiMedia'):
        if mode == 'multiMedia':
            if not multiDados:
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

            else:
                temp = 0
                Vars = []
                SQR = 0
                gl1 = len(multiDados) - 1
                glr = -len(multiDados)
                for k in range(len(multiDados)):
                    Vars.append(np.var(multiDados[k], ddof=1))
                    SQR += (len(multiDados[k]) - 1) * Vars[k]
                    glr += len(multiDados[k])
                    temp += np.mean(multiDados[k]) * len(multiDados[k])
                glt = gl1 + glr
                X = temp / (glr + len(multiDados))
                SQT = 0
                for k in range(len(multiDados)):
                    for n in range(len(multiDados[k])):
                        SQT += (multiDados[k][n] - X)**2
                SQ1 = SQT - SQR
                QM1 = SQ1 / gl1
                QMR = SQR / glr
                QMT = SQT / glt

            

            # else:
            #     gl1 = len(multiDados)-1
            #     glr = len(multiDados) * (len(multiDados[0])-1)
            #     glt = len(multiDados[0]) * len(multiDados) - 1

            #     SQR = 0
            #     SQT = 0
            #     temp = 0
            #     for i in range(len(multiDados)):
            #         media = np.mean(multiDados[i])
            #         for j in range(len(multiDados[i])):
            #             SQR += (multiDados[i][j] - media)**2
            #             temp += multiDados[i][j]
            #     X = temp / (glt + 1)
            #     for i in range(len(multiDados)):
            #         media = np.mean(multiDados[i])
            #         for j in range(len(multiDados[i])):
            #             SQT += (multiDados[i][j] - X)**2
            #     SQ1 = SQT - SQR
            #     QM1 = SQ1 / gl1
            #     QMR = SQR / glr
            #     QMT = SQT / glt

            F1crit = f.ppf(1 - self.alfa, gl1, glr)
            F1calc = QM1 / QMR
            
            pVal = f.cdf(F1calc, gl1, glr)


            print('| Fonte        SQ      gl       QM    |')
            print('| Entre      {:.2f}     {:.0f}       {:.3f} |'.format(SQ1, gl1, QM1))
            print('| Residual   {:.2f}     {:.0f}      {:.3f} |'.format(SQR, glr, QMR))
            print('| Total      {:.2f}    {:.0f}      {:.3f} |'.format(SQT, glt, QMT))
            print()
            print("Fcalc = {:.3f}".format(F1calc))
            print("Fcrit = {:.3f}".format(F1crit))
            print("P-Valor = {:.3f}".format(pVal))
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

        elif mode == 'doisFatRep':
            gl1 = n1 - 1
            gl2 = n2 - 1
            gli = gl1 * gl2
            glr = n1 * n2 * (repeticao - 1)
            glt = n1 * n2 * repeticao - 1

            if QM1:
                SQ1 = QM1 * gl1
            if QM2:
                SQ2 = QM2 * gl2
            if QMI:
                SQI = QMI * gli
            if QMR:
                SQR = QMR * glr
            if QMT:
                SQT = QMT * glt
            if not SQ1:
                SQ1 = SQT - SQR - SQ2 - SQI
            elif not SQ2:
                SQ2 = SQT - SQR - SQ1 - SQI
            elif not SQI:
                SQI = SQT - SQR - SQ1 - SQ2
            elif not SQR:
                SQR = SQT - SQ1 - SQ2 - SQI
            elif not SQT:
                SQT = SQ1 + SQ2 + SQR - SQI

            QM1 = SQ1 / gl1
            QM2 = SQ2 / gl2
            QMI = SQI / gli
            QMR = SQR / glr
            QMT = SQT / glt

            F1calc = QM1 / QMR
            F2calc = QM2 / QMR
            F12calc = QMI / QMR

            F1crit = f.ppf(1 - self.alfa, gl1, glr)
            F2crit = f.ppf(1 - self.alfa, gl2, glr)
            F12crit = f.ppf(1 - self.alfa, gli, glr)

            pVal1 = 1 - f.cdf(F1calc, gl1, glr)
            pVal2 = 1 - f.cdf(F2calc, gl2, glr)
            pVal12 = 1 - f.cdf(F12calc, gli, glr)

            print('| Fonte       SQ          gl        QM          F        Fcrit      p-Valor|')
            print('| Dado 1      {:.2f}       {:.0f}       {:.3f}  {:.3f}  {:.3f}     {:.4f} |'.format(SQ1, gl1, QM1, F1calc, F1crit, pVal1))
            print('| Dado 2      {:.2f}        {:.0f}       {:.3f}   {:.3f}   {:.3f}     {:.4f} |'.format(SQ2, gl2, QM2, F2calc, F2crit, pVal2))
            print('| Interação   {:.2f}         {:.0f}       {:.3f}     {:.3f}    {:.3f}     {:.4f} |'.format(SQI, gli, QMI, F12calc, F12crit, pVal12))
            print('| Residual    {:.2f}         {:.0f}      {:.3f} |'.format(SQR, glr, QMR))
            print('| Total       {:.2f}       {:.0f}      {:.3f} |'.format(SQT, glt, QMT))
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

    def TC(self, SQ1 = None, SQR = None, SQT = None, QM1 = None, QMR = None, QMT = None, n = None):
        somaX = 0
        somaY = 0
        somaXY = 0
        somaX2 = 0
        somaY2 = 0
        if SQ1 or SQR or SQT or QM1 or QMR or QMT:
            gl1 = 1
            glr = n-2
            glt = n-1
            if not SQ1:
                SQ1 = SQT - SQR
            if not SQR:
                SQR = SQT - SQ1
            if not SQT:
                SQT = SQ1 + SQR
            if QM1:
                SQ1 = QM1 * gl1
            if QMR:
                SQR = QMR * glr
            if QMT:
                SQT = QMT * glt
            QM1 = SQ1 / gl1
            QMR = SQR / glr
            QMT = SQT / glt
            
            R2 = SQ1 / SQT
            r = np.sqrt(R2)
            print('| Fonte        SQ        gl          QM    |')
            print('| Entre     {:.2f}     {:.0f}       {:.3f} |'.format(SQ1, gl1, QM1))
            print('| Residual   {:.2f}    {:.0f}      {:.3f} |'.format(SQR, glr, QMR))
            print('| Total      {:.2f}    {:.0f}      {:.3f} |'.format(SQT, glt, QMT))
            print()
            print('|r| = {:.3f}'.format(r))
        else:
            for i in range(len(self.A)):
                somaX += self.A[i]
                somaY += self.B[i]
                somaXY += self.A[i] * self.B[i]
                somaX2 += self.A[i]**2
                somaY2 += self.B[i]**2
            SXX = somaX2 - somaX**2 / self.Na
            SYY = somaY2 - somaY**2 / self.Nb
            SXY = somaXY - somaX * somaY / self.Na

            r = SXY / np.sqrt(SXX * SYY)
            R2 = r**2

            print('r = {:.3f}'.format(r))
        print('R² = {:.3f}'.format(R2))

        tcalc = r * np.sqrt((n - 2) / (1 - R2))
        tcrit = t.ppf(1-self.alfa/2, n- 2)

        if abs(tcalc) > tcrit:
            print("Rejeito H0, há correlação linear")
            return False
        else:
            print("Aceito H0, não há correlação linear")
            return True

    # def TC(self):
    #     somaX = 0
    #     somaY = 0
    #     somaXY = 0
    #     somaX2 = 0
    #     somaY2 = 0
    #     for i in range(len(self.A)):
    #         somaX += self.A[i]
    #         somaY += self.B[i]
    #         somaXY += self.A[i] * self.B[i]
    #         somaX2 += self.A[i]**2
    #         somaY2 += self.B[i]**2
    #     SXX = somaX2 - somaX**2 / self.Na
    #     SYY = somaY2 - somaY**2 / self.Nb
    #     SXY = somaXY - somaX * somaY / self.Na

    #     r = SXY / np.sqrt(SXX * SYY)
    #     R2 = r**2

    #     print('r = {:.3f}'.format(r))
    #     print('R² = {:.3f}'.format(R2))

    #     tcalc = r * np.sqrt((self.Na - 2) / (1 - R2))
    #     tcrit = t.ppf(1-self.alfa/2, self.Na - 2)

    #     if abs(tcalc) > tcrit:
    #         print("Rejeito H0, há correlação linear")
    #         return False
    #     else:
    #         print("Aceito H0, não há correlação linear")
    #         return True
        
        

    def regressao(self, xEmAnalise, n=None, somaX=None, somaX2=None, somaY=None, somaY2=None, somaXY=None, mode = "previsao", Var = None):
        X = somaX / n
        Y = somaY / n
        Sxx = somaX2 - somaX**2 / n
        Sxy = somaXY - somaX * somaY / n
        Syy = somaY2 - somaY**2 / n

        b1 =  Sxy / Sxx
        b0 = Y - b1 * X
        if b0 > 0: print("Y = {:.3f}x + {:.3f}".format(b1, b0))
        else: print("Y = {:.3f}X  {:.3f}".format(b1, b0))

        if mode == "previsao":
            varParcial = np.sqrt(1 + 1/n + (xEmAnalise - X)**2 / Sxx)
        else:
            varParcial = np.sqrt(1/n + (xEmAnalise - X)**2 / Sxx)

        Sr = (Syy - b1 * Sxy) / (n-2)

        # print(Sr)
        # print(varParcial)
        Tcrit = abs(t.ppf(self.alfa / 2, n-2))
        # print("T = ", Tcrit)
        YnoPonto = b0 + b1 * xEmAnalise
        e0 = Tcrit * Sr * varParcial if not Var else Tcrit * Var
        
        print("IC: {:.3f} ± {:.3f}".format(YnoPonto , e0))
        print("IC: [{:.3f}, {:.3f}]".format(YnoPonto - e0, YnoPonto + e0))

        r = Sxy / np.sqrt(Sxx * Syy)
        print("r = {:.3f}".format(r))
        print("R² = {:.3f}".format(r**2))


        print("Teste de Hipótese para a regressão")
        Tb1calc = r * np.sqrt((n-2) / (1-r**2))
        print("Tcalc = {:.3f}".format(Tb1calc))
        print("Tcrit = {:.3f}".format(Tcrit))
        if -Tcrit < Tb1calc < Tcrit:
            return print("Rejeito H0, não há regressão")
        else:
            return print("Aceito H0, há regressão")