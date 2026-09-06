# Átadás-átvételi és független validálási dokumentum

## MethodsX kézirat és reprodukálhatósági csomag

| Tétel | Adat |
|---|---|
| Projekt | `METHODSX_MOBILITY` |
| Kézirat munkacíme | *Comparing Network Indicators and Graphon Distances for Structural Anomaly Detection: A Reproducible Workflow for Directed Weighted Longitudinal Networks* |
| Célfolyóirat és cikktípus | MethodsX, Method Article |
| Átadó | Zsolt T. Kosztyán |
| Átvevő és független validáló | Kornél Dénes |
| Csomagverzió | 0.2.0 |
| Átadás dátuma | 2026-08-28 |
| Életciklusfázis | benyújtás előtti független validálás és kiadási ellenőrzés |

> **Fontos:** ez a dokumentum egy validálásra átadott kiadási jelöltet ír le. A jelenléte önmagában nem jelenti azt, hogy a független validálás megtörtént. A kézirat csak Kornél dokumentált ellenőrzése és mindkét szerző jóváhagyása után tekinthető benyújtásra késznek.

## 1. Az átadás célja és határa

Az átadás célja annak független ellenőrzése, hogy a csomag:

1. az átadott aggregált éves hálózati adatokból reprodukálja az elemzési eredményeket;
2. helyesen valósítja meg a dokumentált hálózati, graphon-inspirált, anomáliadetekciós és szimulációs eljárásokat;
3. teljesíti a beépített analitikai és technikai elfogadási feltételeket;
4. ugyanazokat a számértékeket közli a táblákban, ábrákban és a kézirat szövegében;
5. hibamentesen előállítja a HTML-kéziratot, a MethodsX DOCX-kéziratot, a grafikus absztraktot és a Cover Lettert.

Az átadás kizárólag a MethodsX-cikkre és a hozzá tartozó reprodukálhatósági csomagra vonatkozik. Nem része az eredeti Applied Network Science-kézirat, a `RESBOOK_EN` könyv, egyéb projektek, a folyóirati beadás végrehajtása vagy bármilyen egyéni szintű adminisztratív adat.

## 2. A tudományos állítás röviden

A workflow két, azonos rangú bizonyítékcsaládot hasonlít össze és integrál:

- az éves hálózati jellemzők közvetlenül értelmezhető változásait; valamint
- a node-aligned, graphon-inspirált kernel-, spektrális és strukturális-szerep távolságokat.

A csatornaspecifikus anomáliapontozás és az ismert igazságú szimulációs kalibráció teszi összehasonlíthatóvá a jeleket. A konvergencia egy szélesebb strukturális anomáliára vonatkozó állítást erősít; az eltérés a változás típusáról, léptékéről vagy specifikációérzékenységéről ad diagnosztikai információt. A graphon-inspirált távolságok nem helyettesítik a hagyományos hálózati mutatókat.

## 3. Átadott anyag és belépési pontok

A validálás elsődleges forrása a GitHubra előkészített ZIP vagy annak kicsomagolt `METHODSX_MOBILITY/` könyvtára. A döntő ellenőrzést friss kicsomagolásból vagy friss Git-klónból kell futtatni, nem a korábbi fejlesztői Dropbox-mappából.

| Útvonal | Szerep |
|---|---|
| `README.md` | projektáttekintés és gyors indítás |
| `HANDOVER.md` | ez az átadási és validálási protokoll |
| `REPRODUCIBILITY.md` | részletes futtatási és elfogadási protokoll |
| `CONTRIBUTING.md` | független ellenőrzés és hibajelentés szabályai |
| `DATA_DICTIONARY.md` | a megosztott adatok és mezők leírása |
| `MethodsX_MOB.Rmd` | a kézirat és az ábrák reprodukálható forrása |
| `run_all.R` | a teljes workflow egyetlen fő belépési pontja |
| `scripts/production/00_config.R` | konfiguráció, seedek és futási módok |
| `scripts/production/01_prepare_empirical.R` | éves hálózatok előkészítése |
| `scripts/production/02_analysis_helpers.R` | közös elemzési függvények |
| `scripts/production/03_empirical_benchmark.R` | empirikus elemzés |
| `scripts/production/04_simulation_benchmark.R` | ismert igazságú szimuláció |
| `scripts/production/05_validate_outputs.R` | 16 analitikai ellenőrzés és provenance |
| `scripts/production/06_graphical_abstract.R` | grafikus absztrakt és export |
| `scripts/production/07_freeze_docx_fields.R` | Word-sablon javítása és mezőaudit |
| `scripts/production/08_cover_letter.R` | reprodukálható Cover Letter |
| `scripts/release_check.R` | önálló kiadásijelölt-ellenőrzés |
| `data/raw/` | megosztható aggregált bemenetek |
| `data/derived/` | átadott publication köztes eredmények |
| `submission/` | Cover Letter szerkeszthető és ellenőrzött változata |

Az `output/` könyvtár a futás során jön létre. A ZIP-ben található `data/derived/` fájlok a gyors kézirat-újrarenderelést segítik, de nem helyettesítik a független, `METHODSX_REUSE_INTERMEDIATE=0` beállítású döntő reprodukciót.

## 4. Adatvédelmi és integritási korlát

A csomag kizárólag aggregált, mikrorégiók közötti éves jelentkezési élsúlyokat és node-metaadatot tartalmazhat. Ne kerüljön a repositoryba vagy hibajegybe:

- pályázói vagy hallgatói szintű rekord;
- közvetlen személyazonosító;
- hozzáférési adat, token vagy lokális Dropbox-metaadat;
- korlátozott adminisztratív állomány.

Ha a validáláshoz ilyen adat látszólag szükséges, a futást meg kell állítani és az átadóval kell egyeztetni. A publikált eredmények reprodukciójának az átadott aggregált adatokból kell sikerülnie.

## 5. Javasolt validálási sorrend

### 5.1. Tiszta munkakönyvtár és környezet rögzítése

1. Csomagold ki a ZIP-et egy új, üres könyvtárba, vagy klónozd frissen a repositoryt.
2. Nyisd meg a `METHODSX_MOB.Rproj` fájlt, vagy állítsd a munkakönyvtárat a repository gyökerére.
3. Rögzítsd az operációs rendszert, az R-verziót és a Pandoc-verziót.
4. Telepítsd a hiányzó függőségeket:

```r
source("scripts/install_dependencies.R")
```

Ajánlott környezet: R 4.4.0 vagy újabb és működő Pandoc. A közvetlen R-függőségeket a `DESCRIPTION` sorolja fel. A platformok közötti, bájtszinten azonos grafika nem követelmény; az analitikai CSV-k, a validálási eredmények és a tudományos következtetések egyezése a döntő.

### 5.2. Gyors diagnosztikai futás

Ez egy füstpróba, nem a végső reprodukció:

```r
Sys.setenv(
  METHODSX_MODE = "diagnostic",
  METHODSX_REUSE_INTERMEDIATE = "0",
  METHODSX_RENDER_RMD = "0",
  METHODSX_SKIP_SIMULATION = "0",
  METHODSX_CUT_RESTARTS = "20"
)

source("run_all.R")
```

Elvárt eredmény: a pipeline analitikai hiba nélkül befejeződik. Ha itt hiba van, a teljes publication futást nem érdemes elindítani a hiba dokumentálása előtt.

### 5.3. Döntő, tiszta publication reprodukció

Indíts új R-munkamenetet, majd futtasd a teljes elemzést köztes eredmény újrahasznosítása nélkül:

```r
Sys.setenv(
  METHODSX_MODE = "publication",
  METHODSX_REUSE_INTERMEDIATE = "0",
  METHODSX_RENDER_RMD = "1",
  METHODSX_SKIP_SIMULATION = "0",
  METHODSX_CUT_RESTARTS = "200",
  METHODSX_EMPIRICAL_CUT_SAMPLES = "1000",
  METHODSX_SIM_REPS = "200",
  METHODSX_SIM_CUT_SAMPLES = "128"
)

source("run_all.R")
source("scripts/release_check.R")
warnings()
```

A döntő futás fix seedje `20260817`. A teljes publication futás számításigényes lehet. A konzol teljes kimenetét és a `warnings()` eredményét mentsd el. Váratlan warningot ne minősíts automatikusan ártalmatlannak; add át a pontos szövegét és a futási környezetet.

### 5.4. Automatizált elfogadási feltételek

A döntő futás csak akkor tekinthető sikeresnek, ha:

- `output/diagnostics/METHODSX_PIPELINE_validation.csv` pontosan 16 sort tartalmaz és minden `passed` érték `TRUE`;
- `source("scripts/release_check.R")` hiba nélkül, `Release checks passed.` üzenettel fejeződik be;
- 16 aktuális ábrafájl jön létre: PDF és PNG a grafikus absztrakthoz, a hat fő ábrához és az S1 kiegészítő ábrához;
- létrejön a HTML-kézirat, a DOCX-kézirat és a Cover Letter;
- a végső DOCX-ben nincs frissíthető Word-mező, külső fájlcél vagy érvénytelen relationship ID;
- az output-manifest nem tartalmazza önmagát.

Kulcskimenetek:

```text
output/MethodsX_MOB_figures.html
output/MethodsX_MOB.docx
output/Cover_Letter_MethodsX.docx
output/MethodsX-reference-clean.docx
output/figures/FIG00_graphical_abstract.pdf
output/figures/FIG00_graphical_abstract.png
output/diagnostics/METHODSX_PIPELINE_validation.csv
output/diagnostics/METHODSX_PIPELINE_configuration.csv
output/diagnostics/METHODSX_PIPELINE_sessionInfo.txt
output/diagnostics/METHODSX_FIGURE_manifest.csv
output/diagnostics/METHODSX_PIPELINE_output_manifest.csv
```

## 6. Rögzített numerikus ellenőrzési pontok

Az alábbi értékek a 0.2.0 publication eredményei. Eltérés esetén először a konfigurációt és a bemeneteket kell ellenőrizni. Ha a tiszta, dokumentált futás eredménye mégis eltér, a kéziratot nem szabad változatlanul jóváhagyni.

| Ellenőrzési pont | Elvárt érték |
|---|---:|
| Évek | 2006–2024, megszakítás nélkül |
| Node-ok száma évente | 175 |
| Egymást követő átmenetek | 18 |
| 2016 teljes élsúlya | 358367 |
| 2019 teljes élsúlya | 384035 |
| Becsülők | `block`, `smooth`, `spectral` |
| Fix felbontások | 2, 4, 8, 13 |
| Maximális kalibrált null tévespozitív arány | 0.05, és mindenképpen ≤ 0.10 |
| Közepes hatású directional shift, raw cut lower-bound detekció | 0.860 |
| Közepes hatású local block shock, raw cut lower-bound detekció | 0.615 |
| Közepes hatású role reassignment, raw cut lower-bound detekció | 0.595 |
| Legerősebb empirikus hatperspektívás konszenzus | 2017–2018 |
| 2017–2018 konszenzuspontszám | 0.848 |
| 2017–2018 perspektívák a 75. percentilis felett vagy azon | 5/6 |
| 2015–2016 helyezése konszenzus szerint | 11. |

A 75. percentilis a kéziratban leíró navigációs határ, nem statisztikai szignifikanciaküszöb. A cut-eljárás mindkét előjelen futó alternating maximizationt, determinisztikus és véletlen újraindításokat, valamint egy legacy random floor elemet használ. Az eredmény reprodukálható alsókorlát-becslés; nem egzakt cut-norma és nem relabeling-optimalizált, címkézetlen graphon cut-távolság.

## 7. Kötelező kód- és módszeraudit

Kornél független ellenőrzése legalább az alábbiakat fedje le:

- a környezeti változók ténylegesen a dokumentált konfigurációt állítják-e be;
- a 19 aggregált éves éllista és a node-regiszter helyesen épül-e 175 node-os, node-aligned hálózatokká;
- a súlyozott és irányított hálózati mutatók definíciója megfelel-e a kéziratnak;
- a graphon-inspirált raw intensity és unit-mean shape összehasonlítás nem keveredik-e össze;
- az alternating cut lower-bound mindkét előjelt, a determinisztikus kezdéseket, a véletlen újraindításokat és a retained random floor elemet megfelelően kezeli-e;
- a szimulációs küszöbök kizárólag a null eloszlásból származnak-e, és helyes-e a tévespozitív arány számítása;
- az estimator- és resolution-robosztussági összegzések visszavezethetők-e az alapértékekre;
- a hat perspektíva aggregációja és egyenlő súlyozása megfelel-e a kéziratnak;
- nincs-e adat- vagy eredményszivárgás a kalibráció és az empirikus értékelés között;
- a seedelés és a fájlírás determinisztikus-e a dokumentált határok között.

Külön figyelmet igényel, hogy a módszer node-aligned longitudinális hálózatokra készült. A node-címkék automatikus felcserélhetőségére vagy általános, unlabeled graphon-távolságra nem szabad következtetni.

## 8. Kézirat- és állításaudit

Válassz ki legalább minden absztraktbeli és konklúzióbeli számállítást, majd vezesd vissza a `data/derived/` CSV-kre. Külön ellenőrizendő:

- a szimulációs detekciós arányok és null tévespozitív arány;
- a 2017–2018-as empirikus konszenzus első helye és 0.848-as értéke;
- az öt a hatból perspektívaállítás;
- a 2015–2016-os átmenet 11. helye;
- a hagyományos mutatók és graphon-inspirált távolságok komplementer, nem helyettesítő értelmezése;
- a korlátok között a node alignment, az aggregált jelentkezési adatok, a küszöb- és aggregációs választások, valamint az alsókorlát-jelleg megfelelő megfogalmazása.

Ha bármely kézirati szám nem vezethető vissza egyértelműen egy generált táblára, az `REVISE`, nem pedig automatikus `GO` döntést jelent.

## 9. Vizuális és Word-technikai ellenőrzés

Nyisd meg a teljes HTML- és DOCX-kéziratot, ne csak az első oldalakat.

### Grafikus absztrakt

- minden doboz, címke és összekötő elem a vásznon belül van;
- a bal alsó kutatási kérdések teljesen látszanak: `WHEN? WHAT?` és `HOW ROBUST?`;
- a hálózat-, idősor-, graphonfelület-, kalibráció- és kontribúciószimbólumok nem takarnak szöveget;
- a piros anomáliajel csak a mini idősordiagramon belül fut;
- az ábra 6.8 hüvelykes Word-megjelenítésnél és kicsinyített nézetben is olvasható;
- a PDF vektoros, a PNG pedig megfelelő felbontású.

### DOCX-kézirat

- a Word javítási figyelmeztetés nélkül nyitja meg;
- nem jelenik meg a „más fájlokra hivatkozó mezők frissítése” kérdés;
- a MethodsX-fejléc minden oldalon megmarad;
- egyetlen ábra, táblázat, képaláírás vagy szövegrész sem lóg ki a margón;
- a grafikus absztrakt nincs számozva, a fő és kiegészítő ábrasorrend következetes;
- a matematikai jelölések, görög betűk, kötőjelek és ékezetek helyesek;
- a DOI- és e-mail-hivatkozások működhetnek, de külső helyi fájlhivatkozás nem maradhat;
- a teljes dokumentum végigolvasva sem tartalmaz hiányzó képet vagy sérült oldaltörést.

### Cover Letter

- létrejön az `output/Cover_Letter_MethodsX.docx`;
- a dátum, cím, szerzők és affiliációk helyesek;
- minden szín és tipográfiai elem hibamentesen renderelődik;
- a benyújtás dátumát, a szerzői jóváhagyást, az exkluzivitást, az eredetiséget és az érdekkonfliktus-nyilatkozatot mindkét szerző külön megerősíti.

## 10. Ismert, már kezelt technikai kockázatok

Ezek nem elfogadott warningok, hanem olyan korábbi hibák, amelyek javítása regressziótesztet igényel:

1. A MethodsX referencia-DOCX fejlécében korábban nem numerikus relationship ID okozott `header1.xml` figyelmeztetéseket. A pipeline most futásidőben készít egy szabványosított `output/MethodsX-reference-clean.docx` másolatot, az eredeti sablont változatlanul hagyva.
2. A Word-ábraszámozás frissíthető `SEQ` mezői korábban általános külsőmező-frissítési kérdést válthattak ki. A pipeline hét automatikus ábraszámmezőt fix szöveggé alakít, majd auditálja a DOCX-et.
3. A DOCX újracsomagolásának meg kell őriznie a `word/document.xml` útvonalat; a ZIP-csomagolás ezért mirror módban történik.
4. A Cover Letter színei érvényes, `#` előtagú hexadecimális értékek.
5. A grafikus absztrakt dekoratív elemei korábban szövegkilógást okozhattak; a felirat olvashatósága minden esetben elsőbbséget élvez az ikonok részletességével szemben.

Ha bármelyik korábbi tünet visszatér, azt regressziós hibaként kell jelenteni.

## 11. Eltérések és hibák jelentése

Egy hibajegy vagy validálási megjegyzés minimális tartalma:

```text
Cím:
Súlyosság: BLOCKER / MAJOR / MINOR / NOTE
Operációs rendszer:
R-verzió:
Pandoc-verzió:
Futtatott parancs és környezeti változók:
Érintett fájl vagy ellenőrzési sor:
Elvárt eredmény:
Tényleges eredmény:
Teljes hiba- vagy warning-szöveg:
Legkisebb reprodukáló példa:
Javasolt javítás vagy nyitott kérdés:
```

Súlyossági értelmezés:

- `BLOCKER`: a tiszta futás, az analitikai validálás, a DOCX megnyitása vagy egy központi állítás reprodukciója sikertelen;
- `MAJOR`: a kimenet létrejön, de módszertani, numerikus vagy érdemi narratív eltérés van;
- `MINOR`: korlátozott vizuális, dokumentációs vagy hordozhatósági hiba;
- `NOTE`: nem blokkoló javaslat vagy tisztázandó megfigyelés.

Ne csatolj korlátozott vagy egyéni szintű adatot. A konzolnapló, a session info, a konfigurációs CSV és a releváns validálási sor általában elegendő.

## 12. Átvételi döntés és lezárás

A validálás végén az alábbi három döntés egyikét kell dokumentálni:

- `GO`: minden automatizált ellenőrzés sikeres, a kulcsszámok reprodukálhatók, nincs érdemi kód- vagy állításeltérés, és a teljes vizuális ellenőrzés megfelelt;
- `REVISE`: a projekt tudományosan életképes, de dokumentált javítás szükséges a benyújtás előtt;
- `BLOCK`: a központi eredmény, a teljes reprodukció vagy a dokumentum integritása nem igazolható.

### Validálói összesítő

| Ellenőrzés | Eredmény | Bizonyíték vagy megjegyzés |
|---|---|---|
| Környezet és függőségek rögzítve | ☐ | |
| Diagnosztikai futás sikeres | ☐ | |
| Tiszta publication futás sikeres | ☐ | |
| Mind a 16 analitikai check sikeres | ☐ | |
| Release check sikeres | ☐ | |
| Kulcsszámok függetlenül visszaellenőrizve | ☐ | |
| Cut lower-bound implementáció auditálva | ☐ | |
| Szimulációs kalibráció auditálva | ☐ | |
| Robosztussági összegzések auditálva | ☐ | |
| HTML és DOCX teljes vizuális ellenőrzése kész | ☐ | |
| Grafikus absztrakt megfelelt | ☐ | |
| Word mező-/külsőhivatkozás-audit megfelelt | ☐ | |
| Cover Letter ellenőrizve | ☐ | |
| CRediT-szerepek és szerzői jóváhagyás egyeztetve | ☐ | |
| Adatvédelmi és megoszthatósági korlát ellenőrizve | ☐ | |

```text
Végső döntés: GO / REVISE / BLOCK

Validáló: Kornél Dénes
Dátum:
Operációs rendszer és R-verzió:
Ellenőrzött commit vagy ZIP SHA-256:
Nyitott hibajegyek:
Összefoglaló megjegyzés:

Aláírás vagy dokumentált elektronikus jóváhagyás:
```

## 13. Benyújtás előtti, validáláson túli kapuk

A technikai `GO` döntés után is szükséges:

- mindkét szerző végleges kézirat-jóváhagyása és felelősségvállalása;
- Kornél tényleges hozzájárulásának megfelelő CRediT-nyilatkozat;
- a Cover Letter nyilatkozatainak ismételt ellenőrzése;
- a kód- és adatlicenc végleges kiválasztása a repository nyilvánossá tétele előtt;
- az alkalmazott AI-eszközök felhasználásának szerzői áttekintése, ellenőrzése és a MethodsX aktuális szabályzata szerinti esetleges közzététele;
- a végső ZIP SHA-256 értékének és a GitHub release/tag azonosítójának rögzítése.

Az AI-eszköz nem szerző. A tudományos állításokért, a kódért, az adatok kezeléséért és a beadott szövegért kizárólag a szerzők vállalnak felelősséget.

