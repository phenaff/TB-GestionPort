library(BondValuation)

ai <- AccrInt(
  StartDate=as.Date("2020-03-30"),
  EndDate=as.Date("2021-03-18"),
  Coup=4.625,
  CpY=1,
  Mat=as.Date("2027-03-30"),
  DCC=3,
  RV=100)

print(ai)

res <- BondVal.Yield(
  CP=126.85,
  SETT=as.Date("2021-03-18"),
  Em=as.Date("2012-03-30"),
  Mat=as.Date("2027-03-30"),
  CpY=1,
  Coup=4.625,
  DCC=3)
  
  print(res)