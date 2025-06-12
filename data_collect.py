import requests
import pandas as pd
import time
import calendar
import logging
from datetime import datetime
from Bio import Entrez
from http.client import IncompleteRead

# -------------------------
# 設定
# -------------------------
EMAIL = "g-kasahara-9b9@eagle.sophia.ac.jp"
USER_AGENT = f"MedBERT-Pretrainer/1.0 (mailto:{EMAIL})"
HEADERS = {"User-Agent": USER_AGENT}
Entrez.email = EMAIL
Entrez.tool = "MedBERT-Pretrainer"
API_SLEEP_TIME = 0.5
PAGE_SIZE = 1000
RETMAX_PUBMED = 500

# 現在日付
TODAY = datetime.now().date()
CURRENT_YEAR = TODAY.year
CURRENT_MONTH = TODAY.month
CURRENT_DAY = TODAY.day

# -------------------------
# ログ設定
# -------------------------
logging.basicConfig(filename="data_collect.log", level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")

# -------------------------
# Europe PMC クエリ（PUB_YEAR を除外）
# -------------------------
QUERY_EPMC = '("infectious disease*" OR "virus*" OR "bacteri*" OR "parasit*" OR "fung*")'

# -------------------------
# PubMed クエリ
# -------------------------
BASE_QUERY_PUBMED = f'("infectious disease*" OR "virus*" OR "bacteri*" OR "parasit*" OR "fung*") AND (2015[PDAT] : {CURRENT_YEAR}[PDAT])'

# -------------------------
# 安全なEntrez.read（IncompleteRead 対策）
# -------------------------
def safe_entrez_read(handle, max_retries=3):
    for attempt in range(max_retries):
        try:
            return Entrez.read(handle)
        except IncompleteRead:
            logging.warning(f"IncompleteRead（再試行 {attempt+1}/{max_retries}）")
            time.sleep(2)
    raise RuntimeError("Entrez.read() に失敗しました")

# -------------------------
# Europe PMC 抽出
# -------------------------
def fetch_europe_pmc(query):
    logging.info("[Europe PMC] 取得開始")
    results = []
    base_url = "https://www.ebi.ac.uk/europepmc/webservices/rest/search"
    cursor = "*"
    page_num = 1
    today_str = TODAY.strftime("%Y-%m-%d")

    while True:
        params = {
            "query": f"{query} AND FIRST_PDATE:[2015-01-01 TO {today_str}]",
            "format": "json",
            "pageSize": PAGE_SIZE,
            "cursorMark": cursor,
            "resultType": "core",
        }
        logging.info(f"[Europe PMC] ページ {page_num} 要求中... (cursor: {cursor})")
        logging.debug(f"[Europe PMC] 実クエリ: {params['query']}")
        try:
            response = requests.get(base_url, params=params, headers=HEADERS, timeout=30)
            response.raise_for_status()
            data = response.json()
        except Exception as e:
            logging.error(f"EuropePMC API エラー: {e}")
            break

        page_results = data.get("resultList", {}).get("result", [])
        logging.info(f"[Europe PMC] ヒット件数（このページ）: {len(page_results)}")
        logging.info(f"[Europe PMC] 現在までの累積抽出件数: {len(results)}")

        for entry in page_results:
            doi = entry.get("doi", "").lower()
            title = entry.get("title", "")
            abstract = entry.get("abstractText", "")
            if not abstract:
                logging.debug(f"[除外] abstract が空: DOI={doi}, Title={title}")
                continue
            results.append({
                "DOI": doi,
                "Title": title,
                "Abstract": abstract,
                "Source": "EuropePMC"
            })

        next_cursor = data.get("nextCursorMark")
        if not next_cursor or next_cursor == cursor or not page_results:
            logging.info("[Europe PMC] 取得完了（またはデータ枯渇）")
            break

        cursor = next_cursor
        page_num += 1
        time.sleep(API_SLEEP_TIME)

    logging.info(f"[Europe PMC] 総取得件数: {len(results)}")
    return results


# -------------------------
# PubMed 半月単位の抽出（9999件制限 + 現在日制限）
# -------------------------
def fetch_pubmed_all_by_half_months(query_base, start_year=2015, end_year=CURRENT_YEAR, retmax=500):
    all_results = []
    for year in range(start_year, end_year + 1):
        for month in range(1, 13):
            if year == CURRENT_YEAR and month > CURRENT_MONTH:
                break

            for half in [1, 2]:
                if half == 1:
                    start_day = 1
                    end_day = 15
                else:
                    start_day = 16
                    end_day = calendar.monthrange(year, month)[1]

                if year == CURRENT_YEAR and month == CURRENT_MONTH:
                    if start_day > CURRENT_DAY:
                        continue
                    if end_day > CURRENT_DAY:
                        end_day = CURRENT_DAY

                start_date = f"{year}/{month:02d}/{start_day:02d}"
                end_date = f"{year}/{month:02d}/{end_day:02d}"
                date_query = f"{start_date}[PDAT] : {end_date}[PDAT]"
                full_query = f"{query_base} AND ({date_query})"

                logging.info(f"📅 {year}年{month:02d}月 {half}期 ({date_query})")

                handle = Entrez.esearch(db="pubmed", term=full_query, retmax=0)
                record = safe_entrez_read(handle)
                total_count = int(record["Count"])
                max_limit = 9999

                if total_count == 0:
                    logging.warning("件数 0 のためスキップ")
                    continue
                elif total_count > max_limit:
                    logging.warning(f"件数 {total_count} → 上限 {max_limit} に制限")
                    total_count = max_limit
                else:
                    logging.info(f"件数: {total_count}")

                for start in range(0, total_count, retmax):
                    logging.info(f"[PubMed] {start + 1}/{total_count} 件目から取得開始")
                    try:
                        handle = Entrez.esearch(db="pubmed", term=full_query, retmax=retmax, retstart=start)
                        search_record = safe_entrez_read(handle)
                        id_list = search_record["IdList"]
                        time.sleep(API_SLEEP_TIME)

                        if not id_list:
                            logging.warning("[PubMed] IDリストが空、スキップ")
                            break

                        ids = ",".join(id_list)
                        handle = Entrez.efetch(db="pubmed", id=ids, retmode="xml")
                        records = safe_entrez_read(handle)
                        time.sleep(API_SLEEP_TIME)

                        for article in records.get('PubmedArticle', []):
                            citation = article['MedlineCitation']
                            article_data = citation.get('Article', {})
                            title = article_data.get('ArticleTitle', '')
                            abstract = article_data.get('Abstract', {}).get('AbstractText', [""])
                            abstract_text = " ".join(abstract) if isinstance(abstract, list) else abstract
                            if not abstract_text:
                                continue
                            elocation = article_data.get('ELocationID', [])
                            doi = ""
                            if isinstance(elocation, list):
                                for loc in elocation:
                                    if loc.attributes.get('EIdType') == 'doi':
                                        doi = str(loc).lower()
                            elif hasattr(elocation, 'attributes') and elocation.attributes.get('EIdType') == 'doi':
                                doi = str(elocation).lower()

                            all_results.append({
                                "DOI": doi,
                                "Title": title,
                                "Abstract": abstract_text,
                                "Source": "PubMed"
                            })

                        logging.info(f"[PubMed] 累積取得件数: {len(all_results)}")

                    except Exception as e:
                        logging.error(f"[PubMed] エラー: {e}")
                        continue

    logging.info(f"[PubMed] 全期間の総取得件数: {len(all_results)}")
    return all_results

# -------------------------
# 実行
# -------------------------
def main():
    epmc_data = fetch_europe_pmc(QUERY_EPMC)
    df_epmc = pd.DataFrame(epmc_data)
    df_epmc.to_csv("infectious_disease_epmc.csv", index=False)

    pubmed_data = fetch_pubmed_all_by_half_months(
        query_base=BASE_QUERY_PUBMED,
        start_year=2015,
        end_year=CURRENT_YEAR,
        retmax=RETMAX_PUBMED
    )
    df_pubmed = pd.DataFrame(pubmed_data)
    df_pubmed.to_csv("infectious_disease_pubmed.csv", index=False)

    merged = pd.concat([df_epmc, df_pubmed], ignore_index=True)
    merged = merged.drop_duplicates(subset="DOI", keep="first")
    merged = merged[merged["DOI"].notnull() & (merged["DOI"] != "")]

    merged.to_csv("infectious_disease_merged.csv", index=False)

    print("\n✅ 保存完了:")
    print(" - infectious_disease_epmc.csv")
    print(" - infectious_disease_pubmed.csv")
    print(" - infectious_disease_merged.csv")
    print(f"📊 DOI重複排除後の論文数: {len(merged)}")
    logging.info("保存完了とCSV出力終了")

if __name__ == "__main__":
    main()
