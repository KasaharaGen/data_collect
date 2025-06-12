import requests
import pandas as pd
import time
import calendar
from http.client import IncompleteRead
from Bio import Entrez

# -------------------------
# 設定
# -------------------------
EMAIL = "g-kasahara-9b9@eagle.sophia.ac.jp"
USER_AGENT = f"MedBERT-Pretrainer/1.0 (mailto:{EMAIL})"
HEADERS = {"User-Agent": USER_AGENT}
Entrez.email = EMAIL
Entrez.tool = "MedBERT-Pretrainer"
API_SLEEP_TIME = 1.0
PAGE_SIZE = 1000
RETMAX_PUBMED = 1000
MAX_PUBMED_RECORDS = 10000000

# -------------------------
# クエリ設定
# -------------------------
QUERY_EPMC = """
("infectious disease*" OR "virus*" OR "bacteri*" OR "parasit*" OR "fung*")
AND PUB_YEAR:[2015 TO 2025]
"""

BASE_QUERY_PUBMED = """
(
  "infectious disease"[Title/Abstract] OR
  "communicable disease"[Title/Abstract] OR
  "virus"[Title/Abstract] OR
  "bacteria"[Title/Abstract] OR
  "parasite"[Title/Abstract] OR
  "fungal infection"[Title/Abstract]
)
"""

# -------------------------
# 安全な Entrez.read (IncompleteRead 対応)
# -------------------------
def safe_entrez_read(handle, max_retries=3):
    for attempt in range(max_retries):
        try:
            return Entrez.read(handle)
        except IncompleteRead:
            print(f"⚠️ IncompleteRead → リトライ中 ({attempt + 1}/{max_retries})")
            time.sleep(2)
        except Exception as e:
            print(f"❌ Entrez.read() エラー: {e}")
            break
    raise RuntimeError("Entrez.read() リトライ上限超過")

# -------------------------
# Europe PMC 論文取得
# -------------------------
def fetch_europe_pmc(query):
    print("[Europe PMC] 取得開始")
    results = []
    base_url = "https://www.ebi.ac.uk/europepmc/webservices/rest/search"
    cursor = "*"
    page_num = 1

    while True:
        params = {
            "query": query,
            "format": "json",
            "pageSize": PAGE_SIZE,
            "cursorMark": cursor,
            "resultType": "lite"
        }
        print(f"[Europe PMC] ページ {page_num} 取得中 (カーソル: {cursor})...")
        response = requests.get(base_url, params=params, headers=HEADERS)
        response.raise_for_status()
        data = response.json()
        page_results = data.get("resultList", {}).get("result", [])
        print(f"[Europe PMC] ヒット件数（このページ）: {len(page_results)}")
        for entry in page_results:
            doi = entry.get("doi", "").lower()
            title = entry.get("title", "")
            abstract = entry.get("abstractText", "")
            results.append({
                "DOI": doi,
                "Title": title,
                "Abstract": abstract,
                "FullText": abstract,
                "Source": "EuropePMC"
            })
        next_cursor = data.get("nextCursorMark")
        if not next_cursor or next_cursor == cursor:
            break
        cursor = next_cursor
        page_num += 1
        time.sleep(API_SLEEP_TIME)

    print(f"[Europe PMC] 総取得件数: {len(results)}")
    return results

# -------------------------
# PubMed 論文取得（最大9999件）
# -------------------------
def fetch_limited_articles(query, retmax):
    handle = Entrez.esearch(db="pubmed", term=query, usehistory="y", retmax=0)
    record = safe_entrez_read(handle)
    total_count = int(record["Count"])
    webenv = record["WebEnv"]
    query_key = record["QueryKey"]
    limited_count = min(total_count, 9999)

    if total_count > 9999:
        print(f"⚠️ 件数 {total_count} → 9999件に制限")

    articles = []
    for start in range(0, limited_count, retmax):
        fetch_handle = Entrez.efetch(
            db="pubmed",
            rettype="xml",
            retmode="xml",
            retstart=start,
            retmax=retmax,
            webenv=webenv,
            query_key=query_key,
            timeout=30
        )
        records = safe_entrez_read(fetch_handle)
        time.sleep(API_SLEEP_TIME)

        for article in records.get("PubmedArticle", []):
            citation = article["MedlineCitation"]
            article_data = citation.get("Article", {})
            title = article_data.get("ArticleTitle", "")
            abstract = article_data.get("Abstract", {}).get("AbstractText", [""])
            abstract_text = " ".join(abstract) if isinstance(abstract, list) else abstract
            elocation = article_data.get("ELocationID", [])
            doi = ""
            if isinstance(elocation, list):
                for loc in elocation:
                    if loc.attributes.get("EIdType") == "doi":
                        doi = str(loc).lower()
            elif hasattr(elocation, 'attributes') and elocation.attributes.get("EIdType") == "doi":
                doi = str(elocation).lower()

            articles.append({
                "DOI": doi,
                "Title": title,
                "Abstract": abstract_text,
                "FullText": abstract_text,
                "Source": "PubMed"
            })

    return articles

# -------------------------
# 半月単位に分割取得
# -------------------------
def fetch_pubmed_all_by_half_months_with_limit(base_query, start_year=2015, end_year=2025, retmax=1000):
    all_articles = []
    for year in range(start_year, end_year + 1):
        for month in range(1, 13):
            start_date_1 = f"{year}/{month:02d}/01"
            end_date_1 = f"{year}/{month:02d}/15"
            query_1 = f"{base_query} AND ({start_date_1}[PDAT] : {end_date_1}[PDAT])"
            print(f"\n📅 {year}年{month:02d}月 前半取得中: {query_1}")
            articles_1 = fetch_limited_articles(query_1, retmax)
            print(f"✅ {len(articles_1)} 件取得")
            all_articles.extend(articles_1)

            last_day = calendar.monthrange(year, month)[1]
            start_date_2 = f"{year}/{month:02d}/16"
            end_date_2 = f"{year}/{month:02d}/{last_day}"
            query_2 = f"{base_query} AND ({start_date_2}[PDAT] : {end_date_2}[PDAT])"
            print(f"\n📅 {year}年{month:02d}月 後半取得中: {query_2}")
            articles_2 = fetch_limited_articles(query_2, retmax)
            print(f"✅ {len(articles_2)} 件取得")
            all_articles.extend(articles_2)

    print(f"\n📦 総取得件数（重複含む）: {len(all_articles)}")
    return all_articles

# -------------------------
# 実行関数
# -------------------------
def main():
    #epmc_data = fetch_europe_pmc(QUERY_EPMC)
    #df_epmc = pd.DataFrame(epmc_data)
    #df_epmc.to_csv("infectious_disease_epmc.csv", index=False)

    pubmed_data = fetch_pubmed_all_by_half_months_with_limit(
        base_query=BASE_QUERY_PUBMED,
        start_year=2015,
        end_year=2025,
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

if __name__ == "__main__":
    main()
