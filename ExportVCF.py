import hail as hl
import sys

input_mt_path=sys.argv[1]
output_bucket=sys.argv[2]
identifier=sys.arv[3]

hl.init(spark_conf={
    # Driver settings
    'spark.driver.memory': '32g',
    'spark.driver.cores': '4',

    # Executor settings
    'spark.executor.memory': '20g',
    'spark.executor.cores': '5',
    'spark.executor.instances': '11',

    # Memory management
    'spark.memory.fraction': '0.6',
    'spark.memory.storageFraction': '0.2',

    # Parallelism and shuffle settings
    'spark.default.parallelism': '320',
    'spark.sql.shuffle.partitions': '320',

    # JVM settings
    'spark.executor.extraJavaOptions': '-XX:+UseG1GC -XX:InitiatingHeapOccupancyPercent=35',
    'spark.driver.extraJavaOptions': '-XX:+UseG1GC -XX:InitiatingHeapOccupancyPercent=35',

    # Network and timeout settings (optional but recommended)
    'spark.network.timeout': '800s',
    'spark.executor.heartbeatInterval': '60s'
})

mt = hl.read_matrix_table(input_mt_path)
hl.export_vcf(mt, f"{output_bucket}/{identifier}.vcf.bgz")
