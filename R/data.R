#' Big Five personality questionnaire
#'
#' A dataset containing answers to a Big Five Personality Test conducted by
#' http://personality-testing.info. These data were collected (c. 2012) through
#' an interactive online personality test. The test was constructed with items 
#' from the International Personality Item Pool. Participants were informed that
#' their responses would be recorded and used for research at the beginning of 
#' the test and asked to confirm their consent at the end of the test. The items
#' were rated on a five point scale where 1=Disagree, 3=Neutral, 5=Agree. 
#' All were presented on one page in the order E1, N1, A1, C1, O1, E2,...
#'
#' This dataset is for demonstration purposes only. Please see
#' http://personality-testing.info/privacypolicy.html and 
#' http://personality-testing.info/about for more information.
#'
#' @format A data frame with 19,719 rows of  57 variables:
#' \describe{
#'   \item{race}{1=Mixed Race, 2=Arctic (Siberian, Eskimo), 
#'               3=Caucasian (European), 4=Caucasian (Indian), 
#'               5=Caucasian (Middle East), 6=Caucasian (North African, Other), 
#'               7=Indigenous Australian, 8=Native American, 
#'               9=North East Asian (Mongol, Tibetan, Korean Japanese, etc), 
#'               10=Pacific (Polynesian, Micronesian, etc), 
#'               11=South East Asian (Chinese, Thai, Malay, Filipino, etc), 
#'               12=West African, Bushmen, Ethiopian, 13=Other}
#'  \item{age}{Entered as text (individuals reporting age < 13 were not recorded)}
#'  \item{engnat}{Response to "is English your native language?". 1=yes, 2=no}
#'  \item{gender}{1=Male, 2=Female, 3=Other}
#'  \item{hand}{"What hand do you use to write with?". 1=Right, 2=Left, 3=Both}
#'  \item{country}{The participant's technical location. ISO country code}
#'  \item{source}{How the participant came to the test. Based on HTTP Referer. 
#'                1=from another page on the test website, 2=from google, 
#'                3=from facebook, 4=from any url with ".edu" in its domain name,
#'                6=other source, or HTTP Referer not provided}
#'  \item{E1}{I am the life of the party.}
#'  \item{E2}{I don't talk a lot.}
#'  \item{E3}{I feel comfortable around people.}
#'  \item{E4}{I keep in the background.}
#'  \item{E5}{I start conversations.}
#'  \item{E6}{I have little to say.}
#'  \item{E7}{I talk to a lot of different people at parties.}
#'  \item{E8}{I don't like to draw attention to myself.}
#'  \item{E9}{I don't mind being the center of attention.}
#'  \item{E10}{I am quiet around strangers.}
#'  \item{N1}{I get stressed out easily.}
#'  \item{N2}{I am relaxed most of the time.}
#'  \item{N3}{I worry about things.}
#'  \item{N4}{I seldom feel blue.}
#'  \item{N5}{I am easily disturbed.}
#'  \item{N6}{I get upset easily.}
#'  \item{N7}{I change my mood a lot.}
#'  \item{N8}{I have frequent mood swings.}
#'  \item{N9}{I get irritated easily.}
#'  \item{N10}{I often feel blue.}
#'  \item{A1}{I feel little concern for others.}
#'  \item{A2}{I am interested in people.}
#'  \item{A3}{I insult people.}
#'  \item{A4}{I sympathize with others' feelings.}
#'  \item{A5}{I am not interested in other people's problems.}
#'  \item{A6}{I have a soft heart.}
#'  \item{A7}{I am not really interested in others.}
#'  \item{A8}{I take time out for others.}
#'  \item{A9}{I feel others' emotions.}
#'  \item{A10}{I make people feel at ease.}
#'  \item{C1}{I am always prepared.}
#'  \item{C2}{I leave my belongings around.}
#'  \item{C3}{I pay attention to details.}
#'  \item{C4}{I make a mess of things.}
#'  \item{C5}{I get chores done right away.}
#'  \item{C6}{I often forget to put things back in their proper place.}
#'  \item{C7}{I like order.}
#'  \item{C8}{I shirk my duties.}
#'  \item{C9}{I follow a schedule.}
#'  \item{C10}{I am exacting in my work.}
#'  \item{O1}{I have a rich vocabulary.}
#'  \item{O2}{I have difficulty understanding abstract ideas.}
#'  \item{O3}{I have a vivid imagination.}
#'  \item{O4}{I am not interested in abstract ideas.}
#'  \item{O5}{I have excellent ideas.}
#'  \item{O6}{I do not have a good imagination.}
#'  \item{O7}{I am quick to understand things.}
#'  \item{O8}{I use difficult words.}
#'  \item{O9}{I spend time reflecting on things.}
#'  \item{O10}{I am full of ideas.}
#' }
#' @source \url{http://personality-testing.info/_rawdata/}
"big5"