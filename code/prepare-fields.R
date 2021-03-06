csv <- "UKBB/ukb41181.csv"
colnames <- scan(text = readLines(csv, n = 1), what = "", sep = ",")
available_fields <- setdiff(unique(sub("(.+)-.+\\..+", "\\1", colnames)), "eid")
length(available_fields)  # 2408

dput(sapply(available_fields, function(field) {
  url <- paste0("https://biobank.ctsu.ox.ac.uk/crystal/field.cgi?id=", field)
  gsubfn::strapply(readLines(url, n = 100),
                   "<td>Description:</td><td>(.+)</td>",
                   simplify = 'c')
}))

fields <- c(
  `46` = "Hand grip strength (left)",
  `47` = "Hand grip strength (right)",
  `48` = "Waist circumference",
  `49` = "Hip circumference",
  `50` = "Standing height",
  `102` = "Pulse rate, automated reading",
  `738` = "Average total household income before tax",
  `1070` = "Time spent watching television (TV)",
  `1080` = "Time spent using computer",
  `1160` = "Sleep duration",
  `1180` = "Morning/evening person (chronotype)",
  `1200` = "Sleeplessness / insomnia",
  `1528` = "Water intake",
  `1558` = "Alcohol intake frequency.",
  `1717` = "Skin colour",
  `1727` = "Ease of skin tanning",
  `1747` = "Hair colour (natural, before greying)",
  `2139` = "Age first had sexual intercourse",
  `2178` = "Overall health rating",
  `2267` = "Use of sun/uv protection",
  `2395` = "Hair/balding pattern",
  `2714` = "Age when periods started (menarche)",
  `2754` = "Age at first live birth",
  `3143` = "Ankle spacing width",
  `3144` = "Heel Broadband ultrasound attenuation, direct entry",
  `3146` = "Speed of sound through heel",
  # `3147` = "Heel quantitative ultrasound index (QUI), direct entry",
  `3148` = "Heel bone mineral density (BMD)",
  `3710` = "Length of menstrual cycle",
  `4079` = "Diastolic blood pressure, automated reading",
  `4080` = "Systolic blood pressure, automated reading",
  `5078` = "logMAR in round (left)",
  `5079` = "logMAR in round (right)",
  # `5084` = "Spherical power (right)",
  # `5085` = "Spherical power (left)",
  `6138` = "Qualifications",

  `20015` = "Sitting height",
  `20022` = "Birth weight",
  `20127` = "Neuroticism score",
  `20150` = "Forced expiratory volume in 1-second (FEV1), Best measure",
  `20151` = "Forced vital capacity (FVC), Best measure",
  `20191` = "Fluid intelligence score",
  `20261` = "avMSE",
  `20458` = "General happiness",
  `20459` = "General happiness with own health",
  `21001` = "Body mass index (BMI)",

  `22671` = "Mean carotid IMT (intima-medial thickness) at 120 degrees ",
  `22674` = "Mean carotid IMT (intima-medial thickness) at 150 degrees ",
  `22677` = "Mean carotid IMT (intima-medial thickness) at 210 degrees ",
  `22680` = "Mean carotid IMT (intima-medial thickness) at 240 degrees ",

  `23099` = "Body fat percentage",
  `23100` = "Whole body fat mass",
  `23101` = "Whole body fat-free mass",
  `23102` = "Whole body water mass",
  `23106` = "Impedance of whole body",

  `12336` = "Ventricular rate",
  `12338` = "P duration",
  `12340` = "QRS duration",
  `22330` = "PQ interval",
  `22331` = "QT interval",
  `22332` = "QTC interval",
  `22333` = "RR interval",
  `22334` = "PP interval",

  `30000` = "White blood cell (leukocyte) count",
  `30010` = "Red blood cell (erythrocyte) count",
  `30020` = "Haemoglobin concentration",
  `30030` = "Haematocrit percentage",
  `30040` = "Mean corpuscular volume",
  `30050` = "Mean corpuscular haemoglobin",
  # `30060` = "Mean corpuscular haemoglobin concentration",
  `30070` = "Red blood cell (erythrocyte) distribution width",
  `30080` = "Platelet count",
  `30090` = "Platelet crit",
  `30100` = "Mean platelet (thrombocyte) volume",
  `30110` = "Platelet distribution width",
  `30120` = "Lymphocyte count",
  `30130` = "Monocyte count",
  `30140` = "Neutrophil count",
  # `30150` = "Eosinophil count",
  # `30160` = "Basophil count",
  # `30170` = "Nucleated red blood cell count",
  `30180` = "Lymphocyte percentage",
  `30190` = "Monocyte percentage",
  `30200` = "Neutrophil percentage",
  `30210` = "Eosinophil percentage",
  `30220` = "Basophil percentage",
  # `30230` = "Nucleated red blood cell percentage",
  # `30240` = "Reticulocyte percentage",
  `30250` = "Reticulocyte count",
  `30260` = "Mean reticulocyte volume",
  `30270` = "Mean sphered cell volume",
  `30280` = "Immature reticulocyte fraction",
  # `30290` = "High light scatter reticulocyte percentage",
  `30300` = "High light scatter reticulocyte count",

  `30500` = "Microalbumin in urine",
  `30510` = "Creatinine (enzymatic) in urine",
  `30520` = "Potassium in urine",
  `30530` = "Sodium in urine",

  `30600` = "Albumin",
  `30610` = "Alkaline phosphatase",
  `30620` = "Alanine aminotransferase",
  `30630` = "Apolipoprotein A",
  `30640` = "Apolipoprotein B",
  `30650` = "Aspartate aminotransferase",
  # `30660` = "Direct bilirubin",
  `30670` = "Urea",
  `30680` = "Calcium",
  `30690` = "Cholesterol",
  `30700` = "Creatinine",
  `30710` = "C-reactive protein",
  `30720` = "Cystatin C",
  `30730` = "Gamma glutamyltransferase",
  `30740` = "Glucose",
  `30750` = "Glycated haemoglobin (HbA1c)",
  `30760` = "HDL cholesterol",
  `30770` = "IGF-1",
  `30780` = "LDL direct",
  `30790` = "Lipoprotein A",
  `30800` = "Oestradiol",
  `30810` = "Phosphate",
  `30820` = "Rheumatoid factor",
  `30830` = "SHBG",
  `30840` = "Total bilirubin",
  `30850` = "Testosterone",
  `30860` = "Total protein",
  `30870` = "Triglycerides",
  `30880` = "Urate",
  `30890` = "Vitamin D"
)

cat(purrr::imap_chr(fields, ~ {
  glue::glue(
    '# {.x} ({.y})',
    'y <- fread2(csv, select = "{.y}-0.0")[ind, 1]',
    'SUMMARY(y)',
    'saveRDS(y, "data/ukbb-quant-pheno/.rds")',
    'saveRDS(log(y), "data/ukbb-quant-pheno/log_.rds")',
    .sep = "\n"
  )
}), sep = "\n\n")


fields_binary <- c(
  `1210` = "Snoring",
  `1220` = "Daytime dozing / sleeping (narcolepsy)",
  `1707` = "Handedness (chirality/laterality)",
  `1747` = "Hair colour (natural, before greying)",
  `2139` = "Age first had sexual intercourse",
  `2149` = "Lifetime number of sexual partners",
  # `2159` = "Ever had same-sex intercourse",
  `2207` = "Wears glasses or contact lenses",
  `2237` = "Plays computer games",
  `2247` = "Hearing difficulty/problems",
  `2296` = "Falls in the last year",
  `2463` = "Fractured/broken bones in last 5 years",
  `3710` = "Length of menstrual cycle",
  `3799` = "Headaches for 3+ months",
  `4598` = "Ever depressed for a whole week",
  `20117` = "Alcohol drinker status",
  `20160` = "Ever smoked",
  `20262` = "Myopia diagnosis",
  `20401` = "Ever addicted to any substance or behaviour",
  # `20445` = "Depression possibly related to childbirth",
  `20447` = "Depression possibly related to stressful or traumatic event",
  `20453` = "Ever taken cannabis",
  # `20480` = "Ever self-harmed",
  `20485` = "Ever contemplated self-harm",
  `20510` = "Recent feelings of depression",
  `20511` = "Recent poor appetite or overeating",
  `20512` = "Recent feelings of foreboding",
  `20513` = "Recent thoughts of suicide or self-harm",
  `20533` = "Trouble falling asleep",
  # `20534` = "Sleeping too much",
  # `20535` = "Waking too early",
  `21064` = "Sensitive stomach",
  `21068` = "Diagnosed with coeliac disease or gluten sensitivity",
  NULL
)

cat(purrr::imap_chr(fields_binary, ~ {
  glue::glue(
    '# {.x} ({.y})',
    'y <- fread2(csv, select = "{.y}-0.0")[ind, 1]',
    'y <- ifelse(y < 0, NA, y > 0)',
    'table(y, exclude = NULL)',
    'saveRDS(y, "data/ukbb-binary-pheno/.rds")',
    .sep = "\n"
  )
}), sep = "\n\n")
